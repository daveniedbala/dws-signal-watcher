# DWS signal watcher

Checks raw NOAA and NWS data every 15 minutes for big weather signals that matter to Western New York,
and pushes an alert to Dave's phone (ntfy) when one crosses its line. Results land in `latest.json`,
which the DWS Board reads for its Alerts strip.

Signals: GEFS polar vortex reversal (10 hPa, 60°N), NWS Buffalo winter/wind watches and warnings,
GEFS AO/NAO, CPC 6–10 and 8–14 day temperature outlooks for Buffalo, and Lake Erie water temperature
against 1995–2025. Thresholds are in `THRESHOLDS` at the top of `watch.py`.

Data: NOAA GEFS on AWS Open Data, api.weather.gov, NOAA CPC, NOAA GLERL. Not an official forecast.
