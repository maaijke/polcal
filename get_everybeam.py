import numpy as np
import everybeam
import casacore.tables as pt
from astropy.coordinates import SkyCoord, ITRS
from astropy.time import Time
from pathlib import Path
import shutil

LIBRARY_DIR = Path(__file__).parent


def prepare_ms(newmsname: str, new_dir: SkyCoord, type="hba", tmpdir="./"):
    """Copy template MS to a mock MS
    Args:
        newmsname (str): output msname
        new_dir (astropy.SkyCoord): new phase center
        type: either hba or lba
        tmpdir: str: where to store the new MS
    """

    input_ms = f"{LIBRARY_DIR}/data/template_{type.lower()}.MS.tgz"
    msname = f"{tmpdir}/template_{type.lower()}.MS"
    newmspath = f"{tmpdir}/{newmsname}"
    shutil.unpack_archive(input_ms, tmpdir)
    pt.taql(
        f"SELECT from {msname} where TIME in (SELECT TIME FROM :: LIMIT 1) giving {newmsname} as plain"
    )
    ra = float(new_dir.ra.deg)
    dec = float(new_dir.dec.deg)
    pt.taql(f"UPDATE {newmsname}::FIELD SET DELAY_DIR=[{ra,dec}]deg")
    pt.taql(f"UPDATE {newmsname}::FIELD SET REFERENCE_DIR=[{ra,dec}]deg")
    pt.taql(f"UPDATE {newmsname}::FIELD SET PHASE_DIR=[{ra,dec}]deg")
    pt.taql(f"UPDATE {newmsname}::FIELD SET LOFAR_TILE_BEAM_DIR=[{ra,dec}]deg")
    shutil.rmtree(Path(tmpdir, msname))
    return newmspath


def init_beam(fname):
    b = everybeam.load_telescope(fname)
    phase_dir = pt.table(fname + "/FIELD").getcol("PHASE_DIR")[0][0]
    return b, np.degrees(phase_dir)


def get_itrf_dirs(direction):
    srcpos = direction
    itrf_dirs = srcpos.transform_to(ITRS)
    return itrf_dirs.cartesian.xyz.value.T


def get_station_response_phasecenter(
    mybeam: everybeam.LOFAR, times: Time, freqs: np.ndarray, station_idx=0, rotate=True
):
    data = np.zeros(times.shape + freqs.shape + (2, 2), dtype=complex)
    mstimes = times.mjd * (3600 * 24.0)

    for itime, time in enumerate(mstimes):
        for ifreq, freq in enumerate(freqs):
            data[itime, ifreq] = mybeam.station_response(
                time=time, freq=freq, station_idx=station_idx, rotate=rotate
            )
    return data


def get_station_response(
    mybeam: everybeam.LOFAR,
    times: Time,
    freqs: np.ndarray,
    station_idx=0,
    rotate=True,
    station0=None,
    direction=None,
):
    data = np.zeros(times.shape + freqs.shape + (2, 2), dtype=complex)
    mstimes = times.mjd * (3600 * 24.0)

    for itime, time in enumerate(mstimes):
        for ifreq, freq in enumerate(freqs):
            data[itime, ifreq] = mybeam.station_response(
                time=time,
                freq=freq,
                station_idx=station_idx,
                rotate=rotate,
                station0=station0,
                direction=direction,
            )
    return data


def get_jones(
    station_id: str,
    source: SkyCoord,
    times: Time,
    freqs: np.ndarray,
    tmpdir: str = "./",
    outfile: str = "eb_jones.txt",
    direction=None,
    refdir=None,
    type='hba',
    rotate=True
):
    fname = f"mock_{station_id}.MS"
    msname = prepare_ms(fname, source, tmpdir=tmpdir, type=type)
    stations = pt.table(msname + "/ANTENNA").getcol("NAME")
    stidx = stations.index(f"{station_id}")
    mybeam, _ = init_beam(msname)
    if direction is None or refdir is None:
        jones = get_station_response_phasecenter(
            mybeam, times, freqs, station_idx=stidx, rotate=rotate
        )
    else:
        alldirs = get_itrf_dirs(direction)
        station0 = get_itrf_dirs(refdir)
        jones = get_station_response(
            mybeam,
            times,
            freqs,
            station_idx=stidx,
            rotate=True,
            direction=alldirs,
            station0=station0,
        )
    shutil.rmtree(Path(tmpdir, fname))
    with open(outfile, "w") as newf:
        for iti, ti in enumerate(times.mjd):
            for ifreq, freq in enumerate(freqs):
                myJ = jones[iti, ifreq]
                print(
                    f"{ti} {freq} {myJ[0,0].real} {myJ[0,0].imag} {myJ[0,1].real} {myJ[0,1].imag} {myJ[1,0].real} {myJ[1,0].imag} {myJ[1,1].real} {myJ[1,1].imag}",
                    file=newf,
                )
    return jones

