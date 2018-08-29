from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase
from astroplan import FixedTarget
from astropy.time import Time
import astropy.units as u
import importlib
import pytz
import matplotlib.pyplot as pl
import radvel
import numpy as np
import pandas as pd
import radvel

from astroquery.vizier import Vizier

# catalog_Martins = Vizier.get_catalogs('J/A+A/576/A134')

from toolbox_exop_database.toolbox_exopeu import ExoPEU
from toolbox_exop_database.toolbox_tepcat import ExoPTEPCat
# from toolbox_exop_database.toolbox_composite import ExoP_Composite


import transit_LT

importlib.reload(transit_LT)

def eccentric_anomaly(f, ecc, positive=True):
    """Compute eccentric anomaly from true anomaly.

    :param f: true anomaly in radians
    :param ecc: eccentricity
    :param positive: if True, the result is in [0, 2.pi], else in [-pi, pi].
    :return ee: eccentric anomaly in radians
    """
    ee = 2 * np.arctan(np.tan(f/2) * np.sqrt((1-ecc)/(1+ecc)))  # eccentric anomaly
    if ee < 0 and positive:
        return 2 * np.pi + ee
    else:
        return ee

def mean_anomaly(f, ecc, positive=True):
    """Compute mean anomaly from true anomaly.

    :param f: true anomaly in radians
    :param ecc: eccentricity
    :param positive: if True, the result is in [0, 2.pi], else in [-pi, pi].
    :return m: mean anomaly in radians
    """
    ee = eccentric_anomaly(f, ecc) # eccentric anomaly
    return ee - ecc * np.sin(ee) # mean anomalie

dbEU = ExoPEU()
dbTEP = ExoPTEPCat(catalog_type="planning")

# tz = pytz.timezone("utc")
tz = pytz.timezone("Chile/Continental")

pl_list = ['HD 209458 b', 'HD 189733 b', '51 peg b', 'HD 75289 b', 'tau Boo b', 'HD 179949 b']
# pl_list = ['51 peg b']
orb_phase_list = ["Transit", "Occultation", "Before Occ", "After Occ"]  # First one has to be transit
#orb_phase_list = ["After Occ",]

t_start_sem = Time("2018-08-31 12:00", format='iso', scale='utc')
t_end_sem = Time("2019-10-01 12:00", format='iso', scale='utc')
# t_start_sem = Time("2018-10-01 12:00", format='iso', scale='utc')
# t_end_sem = Time("2019-09-30 12:00", format='iso', scale='utc')
# t_start_sem = Time("2013-06-08 00:00", format='iso', scale='utc')
# t_end_sem = Time("2013-09-30 23:59", format='iso', scale='utc')
t_start = t_start_sem.jd
t_end = t_end_sem.jd

deg_BE_start = (-45 + 180) * u.deg
deg_BE_end = (-15 + 180) * u.deg
deg_AE_start = (15 + 180) * u.deg
deg_AE_end = (45 + 180) * u.deg
deg_around_conj = 7.5

h_min = 1 * u.h
alt_min = 36.03  # 30 <=> airmass=2, 36.03 <=> airmass=1.7

tsamp_per = 1  # %
twilight = False
plot = False

# ephemeris
df_ephem = pd.read_table("/Users/olivier/Softwares/Specific_Analysis/ESPRESSO_GTO/ephemeris/ephemeris.txt", sep="\s+", index_col="planet_name")

dico_properties = {}
dico_windows = {}
file_name = "all_obs_windows_{}.txt"
save_to_file = True

for plnt in pl_list:
    file_name_plnt = (file_name.format(plnt)).replace(" ", "_")
    dico_properties[plnt] = {}
    dico_windows[plnt] = pd.DataFrame()
    pl_properties_EU = dbEU.get_planet_rows(plnt, usematch=True)
    if plnt in []:  # ['HD 209458 b', 'HD 189733 b']
        pl_properties_TEP = dbTEP.get_planet_rows(plnt[:-2].strip(" ").replace(" ", "_"), usematch=True)
        idx_TC = pl_properties_TEP.index[0]
        idx_EU = pl_properties_EU.index[0]
        # dico_properties[plnt]["pl_per"] = pl_properties_TEP.loc[idx_TC, "Period(day)"] * u.d
        # dico_properties[plnt]["pl_omega"] = pl_properties_EU.loc[idx_EU, "omega"] * u.deg
        # dico_properties[plnt]["t_tr"] = pl_properties_TEP.loc[idx_TC, "T0 (HJD or BJD)"] * u.d
        # dico_properties[plnt]["tr_dur"] = pl_properties_TEP.loc[idx_TC, "length"] * u.d
        # dico_properties[plnt]["pl_ecc"] = pl_properties_EU.loc[idx_EU, "eccentricity"]
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["pl_omega"] = pl_properties_EU.loc[idx_EU, "omega"] * u.deg
        dico_properties[plnt]["t_tr"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["tr_dur"] = pl_properties_TEP.loc[idx_TC, "length"] * u.d
        dico_properties[plnt]["pl_ecc"] = pl_properties_EU.loc[idx_EU, "eccentricity"]
        dico_properties[plnt]["transit"] = True
        dico_properties[plnt]["inferior conjunction"] = True
    elif plnt == "HD 189733 b":
        pl_properties_TEP = dbTEP.get_planet_rows(plnt[:-2].strip(" ").replace(" ", "_"), usematch=True)
        idx_TC = pl_properties_TEP.index[0]
        idx_EU = pl_properties_EU.index[0]
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["t_tr"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["pl_omega"] = 90 * u.deg
        dico_properties[plnt]["tr_dur"] = pl_properties_TEP.loc[idx_TC, "length"] * u.d
        dico_properties[plnt]["pl_ecc"] = 0.0
        dico_properties[plnt]["transit"] = True
        dico_properties[plnt]["inferior conjunction"] = True
    elif plnt == "HD 209458 b":
        pl_properties_TEP = dbTEP.get_planet_rows(plnt[:-2].strip(" ").replace(" ", "_"), usematch=True)
        idx_TC = pl_properties_TEP.index[0]
        idx_EU = pl_properties_EU.index[0]
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["t_tr"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["pl_omega"] = 90 * u.deg
        dico_properties[plnt]["tr_dur"] = pl_properties_TEP.loc[idx_TC, "length"] * u.d
        dico_properties[plnt]["pl_ecc"] = 0.0
        dico_properties[plnt]["transit"] = True
        dico_properties[plnt]["inferior conjunction"] = True
    elif plnt == "51 peg b":
        df_51Peg = pd.read_table("/Users/olivier/Data/exoplanet_database/manual_updates/51Peg_v0.txt", comment="#", header=0, sep="\s+")
        # dico_properties[plnt]["pl_per"] = df_51Peg["pl_per"].squeeze() * u.d
        # dico_properties[plnt]["pl_omega"] = df_51Peg["pl_omega"].squeeze() * u.deg
        # dico_properties[plnt]["t_peri"] = (df_51Peg["t_peri"].squeeze() + 2400000) * u.d
        # dico_properties[plnt]["pl_ecc"] = df_51Peg["pl_ecc"].squeeze()
        # dico_properties[plnt]["inferior conjunction"] = False
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["t_ic"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["pl_omega"] = 90 * u.deg
        dico_properties[plnt]["pl_ecc"] = 0.0
        dico_properties[plnt]["transit"] = False
        dico_properties[plnt]["inferior conjunction"] = True
        dico_properties[plnt]["deg_around_conj"] = deg_around_conj
    elif plnt == 'HD 75289 b':
        df_HD752 = pd.read_table("/Users/olivier/Data/exoplanet_database/manual_updates/HD75289_v0.txt", comment="#", header=0, sep="\s+")
        # dico_properties[plnt]["pl_per"] = df_HD752["pl_per"].squeeze() * u.d
        # dico_properties[plnt]["pl_omega"] = df_HD752["pl_omega"].squeeze() * u.deg
        # dico_properties[plnt]["t_peri"] = (df_HD752["t_peri"].squeeze() + 2400000) * u.d
        # dico_properties[plnt]["pl_ecc"] = df_HD752["pl_ecc"].squeeze()
        # dico_properties[plnt]["inferior conjunction"] = False
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["t_ic"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["pl_omega"] = 90 * u.deg
        dico_properties[plnt]["pl_ecc"] = 0.0
        dico_properties[plnt]["transit"] = False
        dico_properties[plnt]["inferior conjunction"] = True
        dico_properties[plnt]["deg_around_conj"] = deg_around_conj
    elif plnt == 'HD 179949 b':
        df_HD179 = pd.read_table("/Users/olivier/Data/exoplanet_database/manual_updates/HD179949_v0.txt", comment="#", header=0, sep="\s+")
        # dico_properties[plnt]["pl_per"] = df_HD179["pl_per"].squeeze() * u.d
        # dico_properties[plnt]["pl_omega"] = df_HD179["pl_omega"].squeeze() * u.deg
        # dico_properties[plnt]["t_peri"] = (df_HD179["t_peri"].squeeze() + 2400000) * u.d
        # dico_properties[plnt]["pl_ecc"] = df_HD179["pl_ecc"].squeeze()
        # dico_properties[plnt]["inferior conjunction"] = False
        dico_properties[plnt]["pl_per"] = df_ephem.loc[plnt, "Period"] * u.d
        dico_properties[plnt]["t_ic"] = df_ephem.loc[plnt, "tic_JDUTC"] * u.d
        dico_properties[plnt]["pl_omega"] = 90 * u.deg
        dico_properties[plnt]["pl_ecc"] = 0.0
        dico_properties[plnt]["transit"] = False
        dico_properties[plnt]["inferior conjunction"] = True
        dico_properties[plnt]["deg_around_conj"] = deg_around_conj
    elif plnt == 'tau Boo b':
        df_tauBoo = pd.read_table("/Users/olivier/Data/exoplanet_database/manual_updates/TauBootis_v0.txt", comment="#", header=0, sep="\s+")
        dico_properties[plnt]["pl_per"] = df_tauBoo["pl_per"].squeeze() * u.d
        dico_properties[plnt]["pl_omega"] = df_tauBoo["pl_omega"].squeeze() * u.deg
        dico_properties[plnt]["t_ic"] = (df_tauBoo["t_ic"].squeeze() + 2400000) * u.d
        dico_properties[plnt]["pl_ecc"] = df_tauBoo["pl_ecc"].squeeze()
        dico_properties[plnt]["transit"] = False
        dico_properties[plnt]["inferior conjunction"] = True
        dico_properties[plnt]["deg_around_conj"] = deg_around_conj
    else:
        dico_properties[plnt]["pl_per"] = pl_properties_EU["pl_per"].squeeze() * u.d
        dico_properties[plnt]["pl_omega"] = pl_properties_EU["pl_omega"].squeeze() * u.deg
        dico_properties[plnt]["t_peri"] = pl_properties_EU["t_peri"].squeeze() * u.d
        dico_properties[plnt]["pl_ecc"] = pl_properties_EU["pl_ecc"].squeeze()
        dico_properties[plnt]["transit"] = False
        dico_properties[plnt]["inferior conjunction"] = False
        dico_properties[plnt]["deg_around_conj"] = deg_around_conj
    target = FixedTarget.from_name(plnt[:-2])
    dico_properties[plnt]["f_t"] = np.pi / 2 - dico_properties[plnt]["pl_omega"].to(u.rad).value
    dico_properties[plnt]["m_t"] = mean_anomaly(dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"])
    target_msg = "#########\nTarget = {}\nPeriod = {}\nEccentricity = {}\nomega = {}\nTransit true anomaly = {} deg\nTransit mean anomaly = {}".format(plnt, dico_properties[plnt]["pl_per"], dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"], np.rad2deg(dico_properties[plnt]["f_t"]), np.rad2deg(dico_properties[plnt]["m_t"]))
    print("\n"+target_msg)
    if save_to_file:
        with open(file_name_plnt, "x") as f:
            f.write(target_msg + "\n")
    target_msg_tr = "Time of inferior conjunction= {} JD_UTC\nComputed from: {}"
    if dico_properties[plnt]["inferior conjunction"]:  # []
        if dico_properties[plnt]["transit"]:
            computed_from = "Transit time"
            dico_properties[plnt]["t0"] = dico_properties[plnt]["t_tr"].to(u.d).value
            dico_properties[plnt]["t_peri"] = radvel.orbit.timetrans_to_timeperi(dico_properties[plnt]["t_tr"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value) * u.day
        else:
            computed_from = "Inferior conjunction time"
            dico_properties[plnt]["t0"] = dico_properties[plnt]["t_ic"].to(u.d).value
            dico_properties[plnt]["t_peri"] = radvel.orbit.timetrans_to_timeperi(dico_properties[plnt]["t_ic"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value) * u.day
    else:
        computed_from = "Periastron passage time"
        dico_properties[plnt]["t0"] = radvel.orbit.timeperi_to_timetrans(dico_properties[plnt]["t_peri"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value) * u.day
    if save_to_file:
        with open((file_name.format(plnt)).replace(" ", "_"), "a") as f:
            f.write(target_msg_tr.format(dico_properties[plnt]["t0"],computed_from) + "\n")
    for key in orb_phase_list:
        dico_properties[plnt][key] = {}
        orbphase_msg = "\n####\nValid observing slots for {}:".format(key)
        print(orbphase_msg)
        if key == "Transit":
            # Transit:
            if dico_properties[plnt]["inferior conjunction"]:  # []
                if dico_properties[plnt]["transit"]:
                    computed_from = "Transit time"
                    orbphase_msg_update = "transit duration = {} h".format(dico_properties[plnt]["tr_dur"].to(u.h).value)
                    # dico_properties[plnt]["t0"] = dico_properties[plnt]["t_tr"].to(u.d).value
                    # dico_properties[plnt]["t_peri"] = radvel.orbit.timetrans_to_timeperi(dico_properties[plnt]["t_tr"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value)
                    dico_properties[plnt][key]["OrbPhase_dur"] = dico_properties[plnt]["tr_dur"]
                    dico_properties[plnt][key]["delta_phi_end"] = (dico_properties[plnt]["tr_dur"] / 2 / dico_properties[plnt]["pl_per"]).decompose().value
                    dico_properties[plnt][key]["delta_phi_start"] = -dico_properties[plnt][key]["delta_phi_end"]
                else:
                    computed_from = "Inferior conjunction time"
                    orbphase_msg_update = "inferior conjunction duration = {} h, corresponding to {} deg".format(2 * dico_properties[plnt]["deg_around_conj"] / 360 * dico_properties[plnt]["pl_per"].to(u.h).value, 2 * dico_properties[plnt]["deg_around_conj"])
                    # dico_properties[plnt]["t0"] = dico_properties[plnt]["t_ic"].to(u.d).value
                    # dico_properties[plnt]["t_peri"] = radvel.orbit.timetrans_to_timeperi(dico_properties[plnt]["t_ic"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value)
                    dico_properties[plnt][key]["OrbPhase_dur"] = (2 * dico_properties[plnt]["deg_around_conj"] *  dico_properties[plnt]["pl_per"] / 360.)
                    dico_properties[plnt][key]["delta_phi_end"] = (dico_properties[plnt]["deg_around_conj"] / 360)
                    dico_properties[plnt][key]["delta_phi_start"] = -dico_properties[plnt][key]["delta_phi_end"]
                dico_properties[plnt][key]["phi_mid"] = 0.
            else:
                computed_from = "Periastron passage time"
                orbphase_msg_update = "inferior conjunction duration = {} h, corresponding to {} deg".format(2 * dico_properties[plnt]["deg_around_conj"] / 360 * dico_properties[plnt]["pl_per"].to(u.h).value, 2 * dico_properties[plnt]["deg_around_conj"])
                # dico_properties[plnt]["t0"] = radvel.orbit.timeperi_to_timetrans(dico_properties[plnt]["t_peri"].to(u.d).value, dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value)
                dico_properties[plnt][key]["OrbPhase_dur"] = (2 * dico_properties[plnt]["deg_around_conj"] *  dico_properties[plnt]["pl_per"] / 360.)
                dico_properties[plnt][key]["phi_mid"] = 0.
                dico_properties[plnt][key]["delta_phi_end"] = (dico_properties[plnt]["deg_around_conj"] / 360)
                dico_properties[plnt][key]["delta_phi_start"] = -dico_properties[plnt][key]["delta_phi_end"]
            dico_properties[plnt][key]["phi_start"] = (dico_properties[plnt][key]["phi_mid"] + dico_properties[plnt][key]["delta_phi_start"]) % 1
            dico_properties[plnt][key]["phi_end"] = (dico_properties[plnt][key]["phi_mid"] + dico_properties[plnt][key]["delta_phi_end"]) % 1
        elif key == "Occultation":
            if dico_properties[plnt]["transit"]:
                orbphase_msg_update = "occultation duration = {} h".format(dico_properties[plnt]["tr_dur"].to(u.h).value)
            else:
                orbphase_msg_update = "superior conjunction duration = {} h, corresponding to {} deg".format(2 * dico_properties[plnt]["deg_around_conj"] / 360 * dico_properties[plnt]["pl_per"].to(u.h).value, 2 * dico_properties[plnt]["deg_around_conj"])
            dico_properties[plnt][key]["phi_mid"] = (mean_anomaly(np.pi + dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"]) - dico_properties[plnt]["m_t"])/ (2 * np.pi)
            orbphase_msg_update += "\nmean anomaly in eclipse = {} in phase".format(dico_properties[plnt][key]["phi_mid"])
            dico_properties[plnt][key]["OrbPhase_dur"] = dico_properties[plnt]["Transit"]["OrbPhase_dur"]
            dico_properties[plnt][key]["delta_phi_end"] = dico_properties[plnt]["Transit"]["delta_phi_end"]
            dico_properties[plnt][key]["delta_phi_start"] = dico_properties[plnt]["Transit"]["delta_phi_start"]
            dico_properties[plnt][key]["phi_start"] = dico_properties[plnt][key]["phi_mid"] + dico_properties[plnt][key]["delta_phi_start"]
            dico_properties[plnt][key]["phi_end"] = dico_properties[plnt][key]["phi_mid"] + dico_properties[plnt][key]["delta_phi_end"]
        elif key == "Before Occ":
            dico_properties[plnt][key]["phi_start"] = ((mean_anomaly(deg_BE_start.to(u.rad).value + dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"]) - dico_properties[plnt]["m_t"]) / (2 * np.pi)) % 1
            dico_properties[plnt][key]["phi_end"] = ((mean_anomaly(deg_BE_end.to(u.rad).value + dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"]) - dico_properties[plnt]["m_t"]) / (2 * np.pi)) % 1
            dico_properties[plnt][key]["OrbPhase_dur"] = (dico_properties[plnt][key]["phi_end"] - dico_properties[plnt][key]["phi_start"]) * dico_properties[plnt]["pl_per"]
            orbphase_msg_update = "phi_start = {}, phi_end = {}\n".format(dico_properties[plnt][key]["phi_start"], dico_properties[plnt][key]["phi_end"])
            orbphase_msg_update += "duration of the period before occultation= {} h".format(dico_properties[plnt][key]["OrbPhase_dur"].to(u.h).value)
        elif key == "After Occ":
            dico_properties[plnt][key]["phi_start"] = ((mean_anomaly(deg_AE_start.to(u.rad).value + dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"]) - dico_properties[plnt]["m_t"])/ (2 * np.pi)) % 1
            dico_properties[plnt][key]["phi_end"] = ((mean_anomaly(deg_AE_end.to(u.rad).value + dico_properties[plnt]["f_t"], dico_properties[plnt]["pl_ecc"]) - dico_properties[plnt]["m_t"])/ (2 * np.pi)) % 1
            dico_properties[plnt][key]["OrbPhase_dur"] = (dico_properties[plnt][key]["phi_end"] - dico_properties[plnt][key]["phi_start"]) * dico_properties[plnt]["pl_per"]
            orbphase_msg_update = "phi_start = {}, phi_end = {}\n".format(dico_properties[plnt][key]["phi_start"], dico_properties[plnt][key]["phi_end"])
            orbphase_msg_update += "duration of the period before occultation= {} h".format(dico_properties[plnt][key]["OrbPhase_dur"].to(u.h).value)
        else:
            raise ValueError
        dico_properties[plnt][key]["cov"] = min(((h_min / dico_properties[plnt][key]["OrbPhase_dur"]).decompose().value, 1.))
        orbphase_msg_update += "\nmin coverage required = {} %\n".format(dico_properties[plnt][key]["cov"] * 100)
        print(orbphase_msg_update)
        if save_to_file:
            with open((file_name.format(plnt)).replace(" ", "_"), "a") as f:
                f.write(orbphase_msg + "\n" + orbphase_msg_update + "\n")
        if plot:
            print("Do plot...")
            fig, ax = pl.subplots()
            ax.set_title("{}: {}".format(plnt, key))
            tt = np.linspace(t_start, t_end, 10000)
            rv_model = radvel.kepler.rv_drive(tt, [dico_properties[plnt]["pl_per"].to(u.d).value, dico_properties[plnt]["t_peri"].to(u.d).value, dico_properties[plnt]["pl_ecc"], dico_properties[plnt]["pl_omega"].to(u.rad).value, 0.5] ) + 0.5
            ax.plot(tt, rv_model)
            # ax.scatter(catalog_Martins[0]["JD"] + 2400000, (catalog_Martins[0]["RV"] - df_51Peg["RV_sys"].squeeze()) / (2 * 1e-3 * df_51Peg["K"].squeeze()) + 0.5)
        else:
            ax = None
        dico_properties[plnt][key]["t_samp"] = tsamp_per/100 * dico_properties[plnt][key]["OrbPhase_dur"].to(u.min).value
        print("Sampling time = {} min".format(dico_properties[plnt][key]["t_samp"]))
        (jds_valid_start,
         jds_valid_end,
         covs) = transit_LT.transit_LT(dico_properties[plnt]["t0"], dico_properties[plnt]["pl_per"].to(u.d).value, target.coord, dico_properties[plnt][key]["phi_start"], dico_properties[plnt][key]["phi_end"], t_start, t_end, observatory="esoparanal", lon=None, lat=None, alt=None,
                                       alt_min=alt_min, tsamp=dico_properties[plnt][key]["t_samp"], cov_min=dico_properties[plnt][key]["cov"], twilight=twilight,
                                       precess=True, nutate=True, aberration=True, refract=True,
                                       verbose=0, plot=plot, ax=ax)
        isotime_valid_start = [Time(x, format= "jd", scale="utc").to_datetime(timezone=tz).strftime("%Y-%m-%d %H:%M") for x in jds_valid_start]
        isotime_valid_end = [Time(x, format= "jd", scale="utc").to_datetime(timezone=tz).strftime("%Y-%m-%d %H:%M") for x in jds_valid_end]
        # isotime_valid_start = [Time(x, format= "jd", scale="utc").jd for x in jds_valid_start]
        # isotime_valid_end = [Time(x, format= "jd", scale="utc").jd for x in jds_valid_end]
        dico_windows[plnt] = pd.concat([dico_windows[plnt], pd.DataFrame({"Type": [key for cov in covs],
                                                                          "OrbPhase dur [h]": [dico_properties[plnt][key]["OrbPhase_dur"] for cov in covs],
                                                                          "cov [%]": [cov * 100 for cov in covs],
                                                                          "ObsWin dur [h]": [cov * dico_properties[plnt][key]["OrbPhase_dur"] for cov in covs],
                                                                          "t_start": jds_valid_start,
                                                                          "t_end": jds_valid_end})],
                                       ignore_index=True)
        if save_to_file:
            print("Writing windows to file {}...".format(file_name_plnt))
            with open(file_name_plnt, "a") as f:
                for cov, start, end in zip(covs, isotime_valid_start,isotime_valid_end):
                    f.write("cov = {:3.0f} %, dur = {:.1f}, {} -> {} [{}]\n".format(cov * 100, cov * dico_properties[plnt][key]["OrbPhase_dur"].to(u.h), start, end, tz))
        else:
            for cov, start, end in zip(covs, isotime_valid_start,isotime_valid_end):
                print("cov = {:3.0f} %, {} -> {} [{}]".format(cov * 100, start, end, tz))
