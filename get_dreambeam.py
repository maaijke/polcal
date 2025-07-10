from dreambeam.rime.scenarios import on_pointing_axis_tracking
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import datetime


def get_jones(
    station_id: str,
    source: SkyCoord,
    times: Time,
    freqs: np.ndarray,
    outfile: str = "db_jones.txt",
    type='hba',
    rotate=True
):
    pointing = [source.ra.rad, source.dec.rad, "J2000"]
    stid = station_id.strip("LBA").strip("HBA")
    timestep = datetime.timedelta(seconds=(times[1] - times[0]).to(u.s).value)
    duration = datetime.timedelta(seconds=(times[-1] - times[0]).to(u.s).value)

    timespy, freqs_db, jones, res = on_pointing_axis_tracking(
        "LOFAR",
        stid,
        type.upper(),
        "Hamaker",
        datetime.datetime.fromisoformat(times[0].isot + "Z"),
        duration,
        timestep,
        pointing,
        do_parallactic_rot=rotate
    )
    freqs_db = np.array(freqs_db)
    freq_index = np.argmin(np.abs(freqs[np.newaxis] - freqs_db[:, np.newaxis]), axis=0)

    with open(f"{outfile}", "w") as myf:
        for iti, ti in enumerate(times.mjd):
            for freq, idx in zip(freqs, freq_index):
                myJ = jones[idx, iti]
                print(
                    f"{ti} {freq} {myJ[0,0].real} {myJ[0,0].imag} {myJ[0,1].real} {myJ[0,1].imag} {myJ[1,0].real} {myJ[1,0].imag} {myJ[1,1].real} {myJ[1,1].imag}",
                    file=myf,
                )
    return jones[freq_index]