#!/usr/bin/env python3
"""DWS signal watcher: checks raw NOAA/NWS data for big weather signals and pings Dave.

Runs every 15 minutes on GitHub Actions (.github/workflows/watch.yml); also runs by hand on the Mac.
Each check has its own interval, so the heavy polar vortex download only runs once a day.
Results go to latest.json (read by the DWS Board build) and alerts.log.

Run by hand:  ~/DWS-Studio/alerts/.venv/bin/python ~/DWS-Studio/alerts/watch.py [--force] [--only vortex]
Thresholds live in THRESHOLDS below. Phone push goes to the ntfy topic in the NTFY_TOPIC env var
(a GitHub secret); config.json can also set ntfy_topic and mac_banner for local runs.
"""
import concurrent.futures as cf
import datetime as dt
import json
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE = HERE / 'state.json'
LATEST = HERE / 'latest.json'
LOG = HERE / 'alerts.log'
UA = 'DWS signal watcher (daveniedbala@gmail.com)'
UTC = dt.timezone.utc

THRESHOLDS = {
    'vortex_steps': [20, 40, 60, 80],   # % of GEFS members with easterly 10 hPa 60N wind; alert on each new step
    'ao_low': -2.0, 'ao_high': 2.5,     # GEFS ensemble-mean AO
    'nao_low': -1.5,                    # GEFS ensemble-mean NAO
    'cpc_strong': 70,                   # CPC probability (%) for Buffalo that counts as strong
    'erie_above_avg': 2.0,
    'noreaster_mb': 990,                # low pressure off the Northeast coast
    'arctic_850c': -10,                 # 850 mb (about 5,000 ft) temperature over Buffalo              # Lake Erie degrees C above the 1995-2025 average for the date
}
INTERVAL_MIN = {'nws': 15, 'cpc': 60, 'aonao': 180, 'erie': 360, 'vortex': 180, 'gefs10': 60}
NWS_EVENTS = re.compile(r'Lake Effect Snow|Winter Storm|Blizzard|Ice Storm|High Wind Warning|Wind Chill Warning|Extreme Cold', re.I)
BUFFALO = (42.89, -78.88)


def get(url, headers=None, timeout=60):
    req = urllib.request.Request(url, headers={'User-Agent': UA, **(headers or {})})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def now():
    return dt.datetime.now(UTC)


def fmt_day(d):
    return d.strftime('%b %-d')


# ---------- checks: each returns (card, [alert messages]) ----------

def check_vortex(state):
    """Share of GEFS members whose zonal-mean 10 hPa wind at 60N turns easterly (a major SSW)."""
    import eccodes
    base = 'https://noaa-gefs-pds.s3.amazonaws.com/gefs.{d}/00/atmos/pgrb2ap5/{m}.t00z.pgrb2a.0p50.f{h:03d}'
    members = ['gec00'] + [f'gep{i:02d}' for i in range(1, 31)]
    leads = list(range(24, 385, 24)) + list(range(432, 841, 48))
    for back in range(3):
        day = (now() - dt.timedelta(days=back)).strftime('%Y%m%d')
        try:
            get(base.format(d=day, m='gep30', h=384) + '.idx', timeout=20)
            break
        except Exception:
            continue
    else:
        raise RuntimeError('no complete GEFS 00z run found')
    run = dt.datetime.strptime(day, '%Y%m%d').replace(tzinfo=UTC)
    if state.get('vortex', {}).get('run') == day and state['vortex'].get('reach') == 35:
        return state['vortex']['card'], []

    def u60(m, h):
        url = base.format(d=day, m=m, h=h)
        try:
            idx = get(url + '.idx', timeout=30).decode().splitlines()
        except Exception:
            if h > 384:  # the 17-35 day part of the run lands a few hours later
                return m, h, None
            raise
        i = next(n for n, l in enumerate(idx) if ':UGRD:10 mb:' in l)
        start = int(idx[i].split(':')[1])
        end = int(idx[i + 1].split(':')[1]) - 1
        g = eccodes.codes_new_from_message(get(url, {'Range': f'bytes={start}-{end}'}, timeout=60))
        try:
            ni = eccodes.codes_get(g, 'Ni')
            lat0 = eccodes.codes_get(g, 'latitudeOfFirstGridPointInDegrees')
            dlat = eccodes.codes_get(g, 'jDirectionIncrementInDegrees')
            row = round(abs(lat0 - 60.0) / dlat)
            vals = eccodes.codes_get_values(g)[row * ni:(row + 1) * ni]
            return m, h, float(vals.mean())
        finally:
            eccodes.codes_release(g)

    winds = {m: {} for m in members}
    with cf.ThreadPoolExecutor(12) as ex:
        for m, h, u in ex.map(lambda a: u60(*a), [(m, h) for m in members for h in leads]):
            if u is not None:
                winds[m][h] = u
    leads = [h for h in leads if all(h in winds[m] for m in members)]
    reach = leads[-1] // 24
    first_rev = {}
    for m, s in winds.items():
        neg = [h for h in leads if s[h] < 0]
        if neg:
            first_rev[m] = neg[0]
    pct = round(100 * len(first_rev) / len(members))
    pct16 = round(100 * sum(1 for h in first_rev.values() if h <= 384) / len(members))
    mean = [round(sum(winds[m][h] for m in members) / len(members), 1) for h in leads]
    earliest = run + dt.timedelta(hours=min(first_rev.values())) if first_rev else None
    hist = state.setdefault('vortex_hist', [])
    hist = [x for x in hist if x['run'] != day][-29:] + [{'run': day, 'pct': pct, 'pct16': pct16}]
    state['vortex_hist'] = hist
    prev = hist[-2]['pct'] if len(hist) > 1 else None
    n = len(first_rev)
    soonest = f', as early as {fmt_day(earliest)}' if earliest else ''
    detail = (f'In {n} of 31 runs of the American weather model, the spinning winds over the Arctic stop and reverse '
              f'within {reach} days{soonest}. That can send Arctic air south 1 to 3 weeks later.'
              + (f' Yesterday: {prev}%.' if prev is not None else ''))
    status = 'alert' if pct >= THRESHOLDS['vortex_steps'][0] else 'watch' if pct >= 10 else 'calm'
    card = {'id': 'vortex', 'icon': '🌪️', 'title': 'Will the polar vortex break down?', 'value': f'{pct}%',
            'status': status, 'detail': detail, 'run': f'GEFS {fmt_day(run)} 00z',
            'series': {'days': [h // 24 for h in leads], 'mean': mean},
            'url': 'https://www.severe-weather.eu/polar-vortex-forecast/'}
    alerts = []
    season = f'{run.year if run.month >= 7 else run.year - 1}'
    sent = state.setdefault('vortex_steps_sent', {}).setdefault(season, [])
    for step in THRESHOLDS['vortex_steps']:
        if pct >= step and step not in sent:
            sent.append(step)
            alerts.append(f'🌪️ The polar vortex may break down. The polar vortex is a huge spinning wind high above the Arctic '
                          f'that usually keeps the coldest air locked up north. In {n} of 31 runs of the American weather model '
                          f'({pct}%), that spin stops and reverses{soonest}. When that happens, Arctic air often spills south '
                          f'1 to 3 weeks later, which can mean cold and snow for Buffalo.'
                          + (f' Yesterday it was {prev}%.' if prev is not None else ''))
    alerts = alerts[-1:]  # one message even if several steps jump at once
    state['vortex'] = {'run': day, 'card': card, 'reach': reach}
    return card, alerts


def check_nws(state):
    data = json.loads(get('https://api.weather.gov/alerts/active?area=NY', {'Accept': 'application/geo+json'}))
    hits = [f['properties'] for f in data['features']
            if 'Buffalo' in (f['properties'].get('senderName') or '') and NWS_EVENTS.search(f['properties']['event'])]
    seen = set(state.get('nws_seen', []))
    alerts = []
    for p in hits:
        if p['id'] in seen:
            continue
        seen.add(p['id'])
        areas = p.get('areaDesc', '')
        meaning = ('A watch means it could get dangerous, so get ready.' if 'Watch' in p['event'] else
                   'A warning means it is happening or about to, so act now.' if 'Warning' in p['event'] else '')
        alerts.append(f"❄️ The National Weather Service in Buffalo just issued a {p['event']} for {areas[:140]}. {meaning}")
    state['nws_seen'] = list(seen)[-200:]
    events = sorted({p['event'] for p in hits})
    card = {'id': 'nws', 'icon': '❄️', 'title': 'Weather Service watches and warnings',
            'value': str(len(hits)) if hits else 'None', 'status': 'alert' if hits else 'calm',
            'detail': ', '.join(events) if hits else 'Nothing for Western New York right now.',
            'url': 'https://www.weather.gov/buf/'}
    return card, alerts


EXPLAIN = {'ao': 'The Arctic Oscillation (AO) measures how tightly cold air is held around the North Pole.',
           'nao': 'The North Atlantic Oscillation (NAO) measures the air pressure pattern over the Atlantic Ocean.'}
TITLES = {'ao': 'Arctic Oscillation: is cold air escaping?', 'nao': 'North Atlantic Oscillation: is cold air steered our way?'}


def check_aonao(state):
    cards, alerts = [], []
    for key, url_part, low, high in [('ao', 'ao.gefs.z1000', THRESHOLDS['ao_low'], THRESHOLDS['ao_high']),
                                     ('nao', 'nao.gefs.z500', THRESHOLDS['nao_low'], None)]:
        rows = get(f'https://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.{url_part}.120days.csv').decode().splitlines()[1:]
        rows = [r.split(',') for r in rows if r]
        init = max(r[2] for r in rows)
        cur = [r for r in rows if r[2] == init]
        by_lead = {}
        for lead, _m, _t, v, valid in cur:
            by_lead.setdefault(int(lead), (valid, []))[1].append(float(v))
        mean = {l: (valid, sum(v) / len(v)) for l, (valid, v) in sorted(by_lead.items())}
        mean = {l: (fmt_day(dt.datetime.strptime(v, '%Y-%m-%d')), x) for l, (v, x) in mean.items()}
        init = fmt_day(dt.datetime.strptime(init, '%Y-%m-%d'))
        lo_l = min(mean, key=lambda l: mean[l][1])
        hi_l = max(mean, key=lambda l: mean[l][1])
        lo_v, hi_v = mean[lo_l][1], mean[hi_l][1]
        name = key.upper()
        status, msg = 'calm', None
        if lo_v <= low:
            status = 'alert'
            msg = (f'🥶 Cold pattern coming? {EXPLAIN[key]} The American model\'s average drops to {lo_v:.1f} by {mean[lo_l][0]}. '
                   + ('A big negative number means the cold air is loosening and can spill south into the eastern U.S.' if key == 'ao' else
                      'A negative number often means a jam in the jet stream that steers cold air into the Great Lakes and East Coast.'))
        elif high is not None and hi_v >= high:
            status = 'alert'
            msg = (f'🌤️ Mild pattern coming? {EXPLAIN[key]} The American model\'s average rises to +{hi_v:.1f} by {mean[hi_l][0]}. '
                   'A big positive number means the cold air stays locked up near the North Pole, so the East tends to stay mild.')
        elif lo_v <= low + 1:
            status = 'watch'
        last = state.get(f'{key}_alerted')
        if msg and last != init:
            prev_state = state.get(f'{key}_status')
            if prev_state != 'alert':
                alerts.append(msg)
            state[f'{key}_alerted'] = init
        state[f'{key}_status'] = status
        cards.append({'id': key, 'icon': '🥶' if lo_v < 0 else '🌤️', 'title': TITLES[key],
                      'value': f'{lo_v:+.1f}', 'status': status,
                      'detail': f'Model average goes as low as {lo_v:+.1f} ({mean[lo_l][0]}) and as high as {hi_v:+.1f}. '
                                + ('Below -2 = cold air can escape south. Above +2.5 = cold stays up north.' if key == 'ao' else
                                   'Below -1.5 = cold air gets steered into the East.'),
                      'series': {'days': list(mean), 'mean': [round(v, 2) for _, v in mean.values()]},
                      'url': 'https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/ao.shtml'})
    return cards, alerts


def check_cpc(state):
    cards, alerts = [], []
    for key, svc, label in [('cpc610', 'cpc_6_10_day_outlk', '6–10 day'), ('cpc814', 'cpc_8_14_day_outlk', '8–14 day')]:
        q = ('https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/' + svc +
             f'/MapServer/0/query?geometry={BUFFALO[1]},{BUFFALO[0]}&geometryType=esriGeometryPoint&inSR=4326'
             '&spatialRel=esriSpatialRelIntersects&outFields=fcst_date,start_date,end_date,prob,cat&returnGeometry=false&f=json')
        feats = json.loads(get(q))['features']
        a = feats[0]['attributes'] if feats else {'cat': 'Normal', 'prob': 33, 'fcst_date': None, 'start_date': None, 'end_date': None}
        cat, prob = a['cat'], round(a['prob'] or 33)
        ms = lambda x: dt.datetime.fromtimestamp(x / 1000, UTC) if x else None
        span = f"{fmt_day(ms(a['start_date']))}–{fmt_day(ms(a['end_date']))}" if a['start_date'] else ''
        issued = a['fcst_date']
        prev = state.get(key)
        word = {'Above': 'warmer than normal', 'Below': 'colder than normal'}.get(cat, 'near normal')
        status = 'alert' if cat in ('Above', 'Below') and prob >= THRESHOLDS['cpc_strong'] else 'watch' if cat == 'Below' else 'calm'
        if prev and issued and prev['issued'] != issued:
            flipped = {prev['cat'], cat} == {'Above', 'Below'}
            strong_new = status == 'alert' and not (prev['cat'] == cat and prev['prob'] >= THRESHOLDS['cpc_strong'])
            if flipped or strong_new:
                was = {'Above': 'warmer than normal', 'Below': 'colder than normal'}.get(prev['cat'], 'near normal')
                alerts.append(f'🗺️ NOAA\'s forecast for {span} {"just flipped" if flipped else "got more confident"}: Buffalo now has a '
                              f'{prob}% chance of being {word} (before: {prev["prob"]}% chance {was}). This comes from NOAA\'s '
                              'Climate Prediction Center, which forecasts whether the weeks ahead will run warmer or colder than usual.')
        if issued:
            state[key] = {'issued': issued, 'cat': cat, 'prob': prob}
        cards.append({'id': key, 'icon': '🌡️', 'title': f'NOAA {label} outlook for Buffalo', 'value': f'{prob}% {cat.lower()}',
                      'status': status, 'detail': f'{span}: {prob}% chance of being {word}.',
                      'url': 'https://www.cpc.ncep.noaa.gov/products/predictions/' + ('610day/' if key == 'cpc610' else '814day/')})
    return cards, alerts


def erie_history():
    cache = HERE / 'history' / 'erie-glsea-1995-2025.json'
    if cache.exists():
        return json.loads(cache.read_text())
    out = {}
    for y in range(1995, 2026):
        listing = get(f'https://apps.glerl.noaa.gov/coastwatch/ftp/glsea/avgtemps/{y}/').decode()
        names = re.findall(r'href="(glsea-temps[^"]+\.dat)"', listing)
        name = next((n for n in names if '1024' in n), names[0])
        for line in get(f'https://apps.glerl.noaa.gov/coastwatch/ftp/glsea/avgtemps/{y}/{name}').decode().splitlines():
            p = line.split()
            if len(p) >= 6 and p[0] == str(y) and p[1].isdigit():
                out.setdefault(p[1], {})[str(y)] = float(p[5])
    cache.parent.mkdir(exist_ok=True)
    cache.write_text(json.dumps(out))
    return out


def check_erie(state):
    lines = get('https://apps.glerl.noaa.gov/coastwatch/ftp/glsea/avgtemps/glsea-temps_1024_3.dat').decode().splitlines()
    rows = [l.split() for l in lines if re.match(r'^\d{4} +\d+', l)]
    y, doy, erie = rows[-1][0], rows[-1][1], float(rows[-1][5])
    hist = erie_history().get(doy, {})
    past = list(hist.values())
    avg = sum(past) / len(past)
    rec_year, rec = max(hist.items(), key=lambda kv: kv[1])
    date = dt.datetime.strptime(f'{y} {doy}', '%Y %j')
    state['erie_c'] = erie
    diff = erie - avg
    f = lambda c: round(c * 9 / 5 + 32)
    record = erie > rec
    status = 'alert' if record or diff >= THRESHOLDS['erie_above_avg'] else 'watch' if diff >= 1 else 'calm'
    alerts = []
    week = date.strftime('%Y-%W')
    if status == 'alert' and state.get('erie_alert_week') != week:
        state['erie_alert_week'] = week
        alerts.append(f'🌊 Lake Erie is {"at a new record" if record else "unusually"} warm for {fmt_day(date)}: '
                      f'{erie * 9 / 5 + 32:.1f}°F, about {diff * 9 / 5:.0f}°F warmer than normal. The record for this date is '
                      f'{rec * 9 / 5 + 32:.1f}°F ({rec_year}). Warm lake water is fuel for lake-effect snow: when cold air blows '
                      'over it, the lake adds heat and moisture, which can build heavier snow bands.')
    card = {'id': 'erie', 'icon': '🌊', 'title': 'Lake Erie water temperature', 'value': f'{f(erie)}°F', 'status': status,
            'detail': f'{fmt_day(date)}: {diff * 9 / 5:.0f}°F warmer than normal ({f(avg)}°F). Record for the date: '
                      f'{rec * 9 / 5 + 32:.1f}°F ({rec_year}). Warmer water = more fuel for lake-effect snow.',
            'url': 'https://coastwatch.glerl.noaa.gov/statistics/average-surface-water-temperature-glsea/'}
    return card, alerts


def event_steps(state, key, pct, steps, reset_below=15):
    """Return the highest new step crossed for an event; forget sent steps once the signal fades."""
    sent = state.setdefault('steps_sent', {}).setdefault(key, [])
    if pct < reset_below:
        sent.clear()
        return None
    new = [s for s in steps if pct >= s and s not in sent]
    sent.extend(new)
    return new[-1] if new else None


def check_storms_cold(state):
    """One pass over the GEFS 10-day forecast: nor'easters, Buffalo's first freeze, Arctic air, lake-effect setups."""
    import eccodes
    import numpy as np
    base = 'https://noaa-gefs-pds.s3.amazonaws.com/gefs.{d}/{c}/atmos/pgrb2ap5/{m}.t{c}z.pgrb2a.0p50.f{h:03d}'
    members = ['gec00'] + [f'gep{i:02d}' for i in range(1, 31)]
    leads = list(range(12, 241, 12))
    for back in range(0, 48, 12):
        t = now() - dt.timedelta(hours=back + 5)
        run = t.replace(hour=0 if t.hour < 12 else 12, minute=0, second=0, microsecond=0)
        day, cyc = run.strftime('%Y%m%d'), run.strftime('%H')
        try:
            get(base.format(d=day, c=cyc, m='gep30', h=240) + '.idx', timeout=20)
            break
        except Exception:
            continue
    else:
        raise RuntimeError('no complete GEFS run found')
    key = f'{day}{cyc}'
    if state.get('gefs10', {}).get('run') == key:
        return state['gefs10']['cards'], []

    def grab(m, h):
        url = base.format(d=day, c=cyc, m=m, h=h)
        idx = get(url + '.idx', timeout=30).decode().splitlines()
        out = {}
        for tag in (':PRMSL:mean sea level:', ':TMP:850 mb:', ':TMIN:2 m above ground:'):
            i = next(n for n, l in enumerate(idx) if tag in l)
            a = int(idx[i].split(':')[1])
            b = int(idx[i + 1].split(':')[1]) - 1 if i + 1 < len(idx) else ''
            g = eccodes.codes_new_from_message(get(url, {'Range': f'bytes={a}-{b}'}, timeout=60))
            try:
                ni, nj = eccodes.codes_get(g, 'Ni'), eccodes.codes_get(g, 'Nj')
                v = eccodes.codes_get_values(g).reshape(nj, ni)  # rows from 90N south, columns from 0E east, 0.5 deg
            finally:
                eccodes.codes_release(g)
            row = lambda lat: round((90 - lat) / 0.5)
            col = lambda lon: round((lon % 360) / 0.5)
            if 'PRMSL' in tag:
                box = v[row(45):row(35) + 1, col(-76):col(-65) + 1]
                out['slp'] = float(box.min()) / 100
            elif 'TMP:850' in tag:
                out['t850'] = float(v[row(43), col(-79)]) - 273.15
            else:
                out['tmin'] = (float(v[row(43), col(-79)]) - 273.15) * 9 / 5 + 32
        return m, h, out

    data = {m: {} for m in members}
    with cf.ThreadPoolExecutor(12) as ex:
        for m, h, o in ex.map(lambda a: grab(*a), [(m, h) for m in members for h in leads]):
            data[m][h] = o
    when = lambda h: run + dt.timedelta(hours=h)
    nmem = len(members)

    def first_hit(test, hours=leads):
        hits = {}
        for m in members:
            hs = [h for h in hours if test(data[m][h])]
            if hs:
                hits[m] = hs[0]
        return hits, round(100 * len(hits) / nmem)

    cards, alerts = [], []
    tag = f'GEFS {fmt_day(run)} {cyc}z'

    # Nor'easter: a low under 990 mb off the Northeast coast (35-45N, 76-65W)
    hits, pct = first_hit(lambda o: o['slp'] < THRESHOLDS['noreaster_mb'])
    deepest = min(data[m][h]['slp'] for m in members for h in leads)
    soon = f'around {fmt_day(when(sorted(hits.values())[len(hits) // 2]))}' if hits else ''
    step = event_steps(state, 'noreaster', pct, [30, 50, 70])
    if step:
        alerts.append(f'🌀 Nor\'easter signal. A nor\'easter is a strong storm that spins up off the East Coast and pulls '
                      f'in wind, rain and snow from the northeast. In {len(hits)} of 31 runs of the American weather model '
                      f'({pct}%), a powerful storm (pressure below {THRESHOLDS["noreaster_mb"]} mb) forms off New England '
                      f'{soon}. Lower pressure means a stronger storm; the strongest run drops to {deepest:.0f} mb.')
    cards.append({'id': 'noreaster', 'icon': '🌀', 'title': "Nor'easter in the next 10 days?", 'value': f'{pct}%',
                  'status': 'alert' if pct >= 30 else 'watch' if pct >= 10 else 'calm',
                  'detail': (f'{len(hits)} of 31 model runs show a strong storm off the Northeast coast {soon}. '
                             if hits else 'No model runs show a strong storm off the Northeast coast. ')
                            + f'Strongest run: {deepest:.0f} mb (lower = stronger).',
                  'run': tag, 'url': 'https://www.wpc.ncep.noaa.gov/'})

    # First freeze at Buffalo (TMIN windows ending at 12z cover the overnight low)
    nights = [h for h in leads if when(h).hour == 12]
    hits, pct = first_hit(lambda o: o['tmin'] <= 32, nights)
    season = str(run.year if run.month >= 7 else run.year - 1)
    frozen = state.get('freeze_seen', {}).get(season)
    if not frozen:
        try:
            obs = json.loads(get(f'https://data.rcc-acis.org/StnData?sid=BUF&sdate={season}-08-01&edate={run:%Y-%m-%d}&elems=mint'))
            frz = [d for d, v in obs['data'] if v not in ('M', 'T', 'S') and float(v.rstrip('A')) <= 32]
            if frz:
                frozen = frz[0]
                state.setdefault('freeze_seen', {})[season] = frozen
        except Exception:
            pass
    first = when(sorted(hits.values())[len(hits) // 2]) if hits else None
    coldest = min(data[m][h]['tmin'] for m in members for h in nights)
    if not frozen and pct >= 50 and state.get('freeze_alerted') != season:
        state['freeze_alerted'] = season
        alerts.append(f'🧊 First freeze coming to Buffalo. Buffalo hasn\'t dropped to 32°F yet this fall. In {len(hits)} of '
                      f'31 runs of the American weather model ({pct}%), it freezes around {fmt_day(first)}. The first freeze '
                      'ends the growing season: frost kills off gardens and tender plants.')
    cards.append({'id': 'freeze', 'icon': '🧊', 'title': 'First freeze in Buffalo', 'value': 'Done' if frozen else f'{pct}%',
                  'status': 'calm' if frozen else 'alert' if pct >= 50 else 'watch' if pct >= 20 else 'calm',
                  'detail': (f'Buffalo first hit 32°F on {fmt_day(dt.datetime.strptime(frozen, "%Y-%m-%d"))}. ' if frozen else
                             f'Not yet this fall. {len(hits)} of 31 model runs freeze within 10 days'
                             + (f', around {fmt_day(first)}. ' if first else '. '))
                            + f'Coldest run: {coldest:.0f}°F.',
                  'run': tag, 'url': 'https://www.weather.gov/buf/'})

    # Arctic air: air about 5,000 ft up colder than -10C (14F) over Buffalo
    hits, pct = first_hit(lambda o: o['t850'] <= THRESHOLDS['arctic_850c'])
    coldest850 = min(data[m][h]['t850'] for m in members for h in leads)
    soon = f'around {fmt_day(when(sorted(hits.values())[len(hits) // 2]))}' if hits else ''
    step = event_steps(state, 'arctic', pct, [50, 80])
    if step:
        alerts.append(f'🥶 Arctic blast heading for Buffalo. In {len(hits)} of 31 runs of the American weather model ({pct}%), '
                      f'air about 5,000 feet up turns colder than 14°F (-10°C) over Buffalo {soon}. Air that cold up high '
                      'usually means a real taste of winter at the ground, and often lake-effect snow.')
    cards.append({'id': 'arctic', 'icon': '🥶', 'title': 'Arctic air over Buffalo?', 'value': f'{pct}%',
                  'status': 'alert' if pct >= 50 else 'watch' if pct >= 20 else 'calm',
                  'detail': f'{len(hits)} of 31 model runs bring air colder than 14°F about 5,000 ft up {soon or "in the next 10 days"}. '
                            f'Coldest run: {coldest850 * 9 / 5 + 32:.0f}°F up there.',
                  'run': tag, 'url': 'https://www.weather.gov/buf/'})

    # Lake-effect setup: lake water at least 13C warmer than the air ~5,000 ft up; snow needs that air below about -5C (23F)
    lake = state.get('erie_c')
    if lake is not None:
        hits, pct = first_hit(lambda o: lake - o['t850'] >= 13 and o['t850'] <= -5)
        _, pct_any = first_hit(lambda o: lake - o['t850'] >= 13)
        soon = f'around {fmt_day(when(sorted(hits.values())[len(hits) // 2]))}' if hits else ''
        step = event_steps(state, 'lakeeffect', pct, [40, 70])
        if step:
            alerts.append(f'❄️ Lake-effect snow setup. When cold air flows over the warmer lake, the lake heats and moistens it, '
                          'and it rises into snow bands. Forecasters look for the lake to be about 23°F (13°C) warmer than the '
                          f'air 5,000 feet up, with that air cold enough for snow. In {len(hits)} of 31 runs of the American '
                          f'weather model ({pct}%), that happens over Lake Erie {soon}. Wind direction decides who gets buried.')
        cards.append({'id': 'lakeeffect', 'icon': '❄️', 'title': 'Lake-effect snow setup?', 'value': f'{pct}%',
                      'status': 'alert' if pct >= 40 else 'watch' if pct >= 15 else 'calm',
                      'detail': (f'{len(hits)} of 31 model runs get cold enough over Lake Erie ({lake * 9 / 5 + 32:.0f}°F water) '
                                 f'for lake-effect snow {soon}. ' if hits else
                                 f'No model runs get cold enough for lake-effect snow yet ({lake * 9 / 5 + 32:.0f}°F water). ')
                                + (f'{pct_any}% show lake-effect rain showers.' if pct_any else ''),
                      'run': tag, 'url': 'https://www.weather.gov/buf/'})

    state['gefs10'] = {'run': key, 'cards': cards}
    return cards, alerts


CHECKS = {'vortex': check_vortex, 'nws': check_nws, 'aonao': check_aonao, 'cpc': check_cpc, 'erie': check_erie, 'gefs10': check_storms_cold}
CARD_IDS = {'vortex': ['vortex'], 'nws': ['nws'], 'aonao': ['ao', 'nao'], 'cpc': ['cpc610', 'cpc814'], 'erie': ['erie'], 'gefs10': ['noreaster', 'freeze', 'arctic', 'lakeeffect']}
ORDER = ['vortex', 'nws', 'aonao', 'cpc', 'erie', 'gefs10']


# ---------- delivery ----------

def notify(messages):
    cfg = json.loads((HERE / 'config.json').read_text()) if (HERE / 'config.json').exists() else {}
    for m in messages:
        with LOG.open('a') as fh:
            fh.write(f'{now().isoformat(timespec="minutes")}  {m}\n')
        if cfg.get('mac_banner') and sys.platform == 'darwin':
            subprocess.run(['osascript', '-e', f'display notification {json.dumps(m)} with title "DWS alert" sound name "Glass"'],
                           check=False)
        topic = os.environ.get('NTFY_TOPIC') or cfg.get('ntfy_topic')
        if topic and not os.environ.get('DWS_QUIET'):
            try:
                urllib.request.urlopen(urllib.request.Request(
                    f"https://ntfy.sh/{topic}", data=m.encode(),
                    headers={'Title': 'DWS alert', 'Tags': 'rotating_light', 'User-Agent': UA}), timeout=30)
            except Exception as e:
                print('ntfy failed:', e, file=sys.stderr)


def main():
    force = '--force' in sys.argv
    only = sys.argv[sys.argv.index('--only') + 1] if '--only' in sys.argv else None
    state = json.loads(STATE.read_text()) if STATE.exists() else {}
    latest = json.loads(LATEST.read_text()) if LATEST.exists() else {'cards': {}, 'recent': []}
    ran = state.setdefault('ran', {})
    new_alerts = []
    for name in ORDER:
        if only and name != only:
            continue
        last = ran.get(name)
        if not force and last and (now() - dt.datetime.fromisoformat(last)).total_seconds() < INTERVAL_MIN[name] * 60:
            continue
        try:
            cards, alerts = CHECKS[name](state)
            for c in cards if isinstance(cards, list) else [cards]:
                c['checked'] = now().isoformat(timespec='minutes')
                latest['cards'][c['id']] = c
            new_alerts += alerts
            ran[name] = now().isoformat()
        except Exception as e:
            print(f'{name} failed: {e!r}', file=sys.stderr)
            for cid, c in latest['cards'].items():
                if cid in CARD_IDS[name]:
                    c['problem'] = f"Couldn't check at {dt.datetime.now().strftime('%-I:%M %p')}"
    for m in new_alerts:
        latest['recent'].insert(0, {'at': now().isoformat(timespec='minutes'), 'text': m})
    latest['recent'] = latest['recent'][:25]
    latest['updated'] = now().isoformat(timespec='minutes')
    STATE.write_text(json.dumps(state, indent=1))
    LATEST.write_text(json.dumps(latest, indent=1))
    if new_alerts:
        notify(new_alerts)
    # Keep the local board's Alerts strip current (skipped on GitHub)
    build = HERE.parent / 'tools' / 'build.mjs'
    if build.exists() and not os.environ.get('CI'):
        subprocess.run(['/usr/local/bin/node', str(build)], check=False, capture_output=True)
    print(f'{len(new_alerts)} new alert(s).', *new_alerts, sep='\n')


if __name__ == '__main__':
    main()
