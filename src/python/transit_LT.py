"""Main module to compute ephemeris.
"""
from numpy import sign, floor, ceil, ones, nan, arange, concatenate, ones_like, logical_and, logical_or, mean, cumsum, array, where
from PyAstronomy.pyasl import eq2hor, sunpos
from astropy.coordinates import get_sun
from astropy.time import Time
import matplotlib.pyplot as pl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import logging
import pytz

tz = pytz.timezone("utc")
# tz = pytz.timezone("Chile/Continental")

# from astropy.coordinates import SkyCoord

logger = logging.getLogger()


def transit_LT(tref, per, coords, phi_start, phi_end, t_start, t_stop, observatory=None, lon=None, lat=None, alt=None,
               alt_min=30., tsamp=15, cov_min=1, twilight=False,
               precess=True, nutate=True, aberration=True, refract=True,
               verbose=0, plot=False, ax=None):
    """Compute the ephemeris and more...

    :param float tref: Reference time  for the computation of the orbital phase in days (reference frame has to be consistent between T0, t_start and t_end)
    :param float per: Orbital period in days
    :param SkyCoord coords: SkyCoord object which contain the coordinates of the target in epoch J2000 FK5
    :param float dec: String giving the declination in degree
    :param float phi_start: Begin of the mean orbital phase range of interest (between 0 and 1). The mean orbital phase is defined using the mean anomaly M, phi = M/(2.Pi). So be careful with eccentric orbits.
    :param float phi_end: End of the mean orbital phase range of interest (between 0 and 1). The mean orbital phase is defined using the mean anomaly M, phi = M/(2.Pi). So be careful with eccentric orbits.
    :param float t_start: Beginning of the observing period in days (reference frame has to be consistent between T0, t_start and t_end)
    :param float t_stop: End of the observing period in days (reference frame has to be consistent between T0, t_start and t_end)
    :param string observatory: A string identifying the observatory. Default is HS for Hamburger Sternwarte. If given, the observer's longitude, latitude, and altitude are set automatically (and must not be given separately then)
    :param float lat: Latitude of the observatory in degrees. Specify South Latitude with a negative sign. Default is the latitude of the Hamburger Sternwarte.
    :param float long: Longitude of the observatory in degrees. Specify West longitude with a negative sign. Default is the longitude of the Hamburger Sternwarte.
    :param float alt: Altitude of teh observatory in meter
    :param float alt_min: minimum local altitude in degrees for the target to be observed in good conditions. Default is 30 which correspond to an airmass of 2.
    :param float tsamp: Time sampling for the computation in minutes, default is 15 min.
    :param float cov_min: Minimum fraction of the orbital phase which has to be covered to be considered valid (between 0 and 1).
    :param bool twilight: If True, twilight is included in the valid observation time.
    :param bool precess: If True (default), Earth's precess motion is considered in the calculations.
    :param bool nutate: If True (default), Earth's nutation is considered in the calculations.
    :param bool aberration: If True (default), the annual aberration is considered in the calculations.
    :param bool refract: If True, the atmospheric refraction is considered in the calculations.
    :param int verbose: If >0, then log the execution
    :param bool plot: If True make a summary plot.
    """
    if verbose:
        logger.debug("Ephemeris of the planet: period={}, t_ref={}".format(per, tref))
        logger.debug("Definition of the full observing period: t_start={}, t_stop={}".format(t_start, t_stop))
    #### Compute the last transit date before the beginning of the observing period
    tref_start = tref + floor((t_start - tref) / per) * per  # Work whatever t_start and t0
    if verbose:
        logger.debug("Last t_ref occurence before the beginning of the observing run: {}".format(tref_start))
    #### Compute number of orbital period until the last transit before the end of the observing period
    nper_obs = int(floor((t_stop - tref_start) / per))
    if nper_obs == 0:
        logger.warning("There is no orbital phase slot ({}, {}) between {} and {}! Return empty lists"
                       "".format(phi_start, phi_end, t_start, t_stop))
        return [], [], []
    if verbose:
        logger.debug("Number of orbital period until the last transit before the end of the observing period: {}".format(nper_obs))
    #### Ensure that the orbital phase range is continuous (phi_end = 1 + phi_end if needed)
    if phi_start > phi_end:
        phi_end += 1  # If phi_start > phi_end it means the the phase range starts before transit and end after, so we convert the phi_end to be above 1.
    if verbose:
        logger.debug("Definition of the orbital phase period of interest: phi_start={}, phi_end={}".format(phi_start, phi_end))
    #### Compute the times corresponding to phi_start and phi_end for all orbital period within the observing period
    t_phirange = ones((nper_obs, 2)) * nan
    t_phirange[:, 0] = tref_start + (phi_start + arange(nper_obs)) * per
    t_phirange[:, 1] = tref_start + (phi_end + arange(nper_obs)) * per
    # print(t_phirange)
    if verbose:
        logger.debug("time of the start and end of each orbital phase of interest during the observing period:\n"
                     "shape of the array: {}\narray:".format(t_phirange.shape, t_phirange))
        for idx_range, t_start_stop in zip(range(nper_obs), t_phirange):
            logger.debug("Phase range {}: {} -> {}".format(idx_range, t_start_stop[0], t_start_stop[1]))
    #### For each phase range get the local horizon coords (alt-az) of the target with a 15 min sampling
    # Compute the Julian dates vector containing the julian dates of the 15 min samplings of all the orbital phase ranges one after the other.
    if verbose:
        logger.debug("Sampling time of the orbital phase periods: {} min".format(tsamp))
    jds = [arange(t_phirange[ii, 0], t_phirange[ii, 1], tsamp / 60. / 24., dtype=float) for ii in range(nper_obs)]
    l_idx_size_range = [int(jd.size) for jd in jds]

    jds = concatenate(jds, axis=0)
    l_idx_range_start = cumsum([0] + l_idx_size_range[:-1])
    l_idx_range_stop = concatenate((cumsum(l_idx_size_range[:-1], dtype=int), [int(jds.size)]), axis=0)
    if verbose:
        logger.debug("Size of the concatenated list of times to evaluate:{}\n".format(len(jds)))
        logger.debug("Description of each orbital phase range:")
        for idx_range, phase_range_size, idx_start, idx_stop in zip(range(nper_obs), l_idx_size_range, l_idx_range_start, l_idx_range_stop):
            logger.debug("Phase range {}: Size = {}; indexes = {} -> {}; times {} -> {}".format(idx_range, phase_range_size, idx_start, idx_stop - 1, jds[idx_start], jds[idx_stop - 1]))
    # Compute the altitude, azimuth and hourangle of the target for each julian date
    alts, azs, has = eq2hor(jds, ones_like(jds) * coords.ra.deg, ones_like(jds) * coords.dec.deg, observatory=observatory, lon=lon, lat=lat, alt=alt, B1950=False)
    if verbose:
        logger.debug("Target Altitude, Azimuth, hourangle computed !")
        logger.debug("Target Altitude: {}".format(alts))
    #### For each phase range get the local horizon coords (alt-az) of the Sun with the same 15 min sampling
    sun_coords = get_sun(Time(jds, format="jd", scale="utc"))
    # ra_sun = []
    # dec_sun = []
    # for ii in range(nper_obs):
    #     _, ra_sun_ii, dec_sun_ii = sunpos(jd=t_phirange[ii, 0], end_jd=t_phirange[ii, 1], jd_steps=int(floor((t_phirange[ii, 1] - t_phirange[ii, 0]) / (tsamp/60./24.))))
    #     ra_sun.append(ra_sun_ii)
    #     dec_sun.append(dec_sun_ii)
    #     print(len(ra_sun_ii), len(ra_sun_ii))
    # ra_sun = concatenate(ra_sun, axis=0)
    # dec_sun = concatenate(dec_sun, axis=0)
    # if verbose:
    #     logger.debug("Size of the concatenated list of RA for the Sun for each time:{}\n"
    #                  "Size of the concatenated list of DEC for the Sun for each time:{}\n"
    #                  "".format(len(ra_sun), len(dec_sun)))
    # alts_sun, azs_sun, has_sun = eq2hor(jds, ra_sun, dec_sun, observatory=observatory, lon=lon, lat=lat, alt=alt, B1950=False)
    alts_sun, azs_sun, has_sun = eq2hor(jds, sun_coords.ra, sun_coords.dec, observatory=observatory, lon=lon, lat=lat, alt=alt, B1950=False)
    if verbose:
        logger.debug("Sun Altitude computed !")
    #### Find periods of : day, twilight, and night
    # day_cond = alts_sun[0] >= 0.
    twi_cond = logical_and(alts_sun > -18., alts_sun < 0.)
    night_cond = alts_sun <= -18.
    #### Find the periods where the target is high enough in the sky
    alt_cond = alts >= alt_min
    #### Find the periods of valid observing time, where the target is high enough and it's dark (night and eventually twilight)
    if twilight:
        darkness_cond = logical_or(twi_cond, night_cond)
    else:
        darkness_cond = night_cond
    valid_cond = logical_and(alt_cond, darkness_cond)
    if verbose:
        logger.debug("Total percentage of the valid observing time = {:4.2f}%".format(valid_cond.mean() * 100))
    #### For each phase range compute the fraction of valid observing time
    #### and the start and end julian dates for valid observation
    coverage = []
    jds_valid_start = []
    jds_valid_end = []
    covs = []
    if verbose:
        logger.debug("Minimal coverage required = {}%".format(cov_min * 100))
    if plot:
        if ax is None:
            fig, ax = pl.subplots()
        ax.set_xlim(t_start, t_stop)
        ax.set_adjustable('box-forced')
        ax_up = ax.twiny()
        ax_up.set_adjustable('box-forced')
        ax_up.set_xlim(t_start, t_stop)
        ax_up.set_xticks(arange(ceil(t_start), floor(t_stop), step=1), minor=True)
        ax_up.set_xticks(arange(ceil(t_start), floor(t_stop), step=10))
        xuptpos = ax_up.get_xticks()
        ax_up.set_xticklabels(map(lambda x: "{}".format(Time(x, format="jd", scale="utc").to_datetime(timezone=tz).strftime("%Y-%m-%d %H:%M")), xuptpos), rotation="vertical")
        ymin = 0
        ymax = 5.5
        ax.set_ylim(ymin, ymax)
        ytpos = [0.5, 2, 3.5, 5]
        ytlabels = ["Orbital phase ranges", "Darkness", "Target Alt", "Total"]
        ax.set_yticks(ytpos)
        ax.set_yticklabels(ytlabels)
        xtpos = ax.get_xticks()
        jd_timeoffset = 2400000
        ax.set_xticklabels(map(lambda x: "{:.0f}".format(x), xtpos - jd_timeoffset), rotation="vertical")
        ax.set_xlabel("time [JD-{:d}]".format(jd_timeoffset))
        pl.vlines(x=t_phirange[:, 0], ymin=ymin, ymax=ymax, color='k', linestyles='dashed', linewidth=0.5)
        pl.vlines(x=t_phirange[:, 1], ymin=ymin, ymax=ymax, color='k', linestyles='dashed', linewidth=0.5)
        boxes_OrbPhase = []
        boxes_darkness = []
        boxes_alt = []
        boxes_valid = []
    for ii in range(nper_obs):
        valid_cond_ii = valid_cond[l_idx_range_start[ii]:l_idx_range_stop[ii]]
        coverage.append(mean(valid_cond_ii))
        if verbose:
            logger.debug("Period {}: Percentage of valid observing time = {:4.2f}%".format(ii, valid_cond_ii.mean() * 100))
        if coverage[ii] >= cov_min:
            covs.append(coverage[ii])
            jds_valids = jds[l_idx_range_start[ii]:l_idx_range_stop[ii]][where(valid_cond_ii)]
            jds_valid_start.append(jds_valids.min())
            jds_valid_end.append(jds_valids.max())
        if plot:
            boxes_OrbPhase.append(Rectangle((t_phirange[ii, 0], 0), t_phirange[ii, 1] - t_phirange[ii, 0], 1))  # Orbital Phases goes from 0 to 1 in y
            for jd, dark, alt, val in zip(jds[l_idx_range_start[ii]:l_idx_range_stop[ii]], darkness_cond[l_idx_range_start[ii]:l_idx_range_stop[ii]], alt_cond[l_idx_range_start[ii]:l_idx_range_stop[ii]], valid_cond_ii):
                if dark:
                    boxes_darkness.append(Rectangle((jd - tsamp / 2 / 60 / 24, 1.5), tsamp / 60 / 24, 1))  # Darkness goes from 1.5 to 2.5 in y
                if alt:
                    boxes_alt.append(Rectangle((jd - tsamp / 2 / 60 / 24, 3), tsamp / 60 / 24, 1))  # Altitude goes from 3 to 4 in y
                if val:
                    boxes_valid.append(Rectangle((jd - tsamp / 2 / 60 / 24, 4.5), tsamp / 60 / 24, 1))  # Total goes from 4.5 to 5.5 in y
    if plot:
        pc_OrbPhase = PatchCollection(boxes_OrbPhase, facecolor="C1", alpha=1, edgecolor="C1")
        pc_Dark = PatchCollection(boxes_darkness, facecolor="C0", alpha=1, edgecolor="C0")
        pc_Alt = PatchCollection(boxes_alt, facecolor="C4", alpha=1, edgecolor="C4")
        pc_Val = PatchCollection(boxes_valid, facecolor="C2", alpha=1, edgecolor="C2")
        ax.add_collection(pc_OrbPhase)
        ax.add_collection(pc_Dark)
        ax.add_collection(pc_Alt)
        ax.plot(jds, (alts / 90) + 3)
        xmin, xmax = ax.get_xlim()
        ax.hlines((30 / 90) + 3, xmin, xmax, linestyles="dashed")
        ax.add_collection(pc_Val)
        pl.tight_layout()
    return jds_valid_start, jds_valid_end, covs
