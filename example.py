from get_dreambeam import get_jones as get_db_jones
from get_everybeam import get_jones as get_eb_jones


from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

freqs = np.arange(100e6, 110e6, 200e6/1024)
startTime = "2025-04-25T03:30:00"
times = Time(startTime) + np.arange(300) * u.s
pointing = [5.069077192739558, 0.38194712622773974]  #b1919 in rad
source = SkyCoord(ra=pointing[0] * u.rad, dec=pointing[1] * u.rad)
stid = "DE601HBA"

db_jones = get_db_jones(station_id=stid, source=source,times=times,freqs=freqs, outfile="test_db.txt") # freq x time x 2 x 2
eb_jones = get_eb_jones(station_id=stid, source=source,times=times,freqs=freqs, outfile="test_eb.txt") # time x freq x 2 x 2