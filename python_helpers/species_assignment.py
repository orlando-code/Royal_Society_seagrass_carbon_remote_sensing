import pandas as pd

# Points data


def assign_points_seagrass_species(seabedtype: str) -> str:
    """Assign seagrass species based on the seabedtype."""
    if seabedtype in [
        "1120 - Posidonia beds (Posidonion oceanicae)",
        "A2.131 - Facies of banks of dead leaves of [Posidonia oceanica] and other phanerogams",
        "A5.535 - [Posidonia] beds",
        'A5.5353 - Facies of dead "mattes" of [Posidonia oceanica] without much epiflora',
        "III.5.1. Association with Posidonia oceanica (P)",
        "III.5.1. Posidonia oceanica meadows (= Association with Posidonia oceanica)",
    ]:
        return "Posidonia oceanica"
    elif seabedtype in [
        "A2.61 - Seagrass beds on littoral sediments",
        "A5.39 - Mediterranean communities of coastal terrigenous muds",
        "A5.53 - Sublittoral seagrass beds",
        "Not Reported",
        "Seagrass beds",
        "unspecified unspecified",
    ]:
        return "Unspecified"
    elif seabedtype in [
        "A2.611 - Mainland Atlantic [Zostera noltii] or [Zostera angustifolia] meadows",
        "A2.6111 - [Zostera noltii] beds in littoral muddy sand",
        "LS.LMp.LSgr.Znol",
        "LS.LMp.LSgr.Znol - Zostera noltii beds in littoral muddy sand",
        "Pontic  Zostera noltii meadows (1-3 m)",
        "A2.6111",
    ]:
        return "Zostera noltei"
    elif seabedtype in [
        "A5.53131 - Association with [Cymodocea nodosa] on well sorted fine sands",
        "A5.53131",
    ]:
        return "Cymodocea nodosa"
    elif seabedtype in [
        "A5.5331 - [Zostera marina]/[angustifolia] beds on lower shore or infralittoral clean or muddy sand",
        "SS.SMp.SSgr.Zmar",
        "SS.SMp.SSgr.Zmar - Zostera marina/angustifolia beds on lower shore or infralittoral clean or muddy sand",
        "Zostera beds",
        "A5.5331",
    ]:
        return "Zostera marina"
    elif seabedtype in [
        "Halophila baillonis",
        "Halophila decipiens",
        "Halophila engelmannii",
        "Halophila spp.",
    ]:
        return "Halophila stipulacea"
    elif seabedtype in [
        "Pontic mixed Zostera noltii- Zannichellia palustris-Zostera marina meadows ( 2-4 m)"
    ]:
        return "Zostera marina and Zostera noltei"
    elif seabedtype in [
        "A5.5343 - [Ruppia maritima] in reduced salinity infralittoral muddy sand",
        "Halodule beaudettei",
        "Halodule wrightii",
        "Syringodium filiforme",
        "Thalassia testudinum",
    ]:
        return "Remove"
    else:
        return "Other"


def process_points_df_species(pts_df) -> pd.DataFrame:
    """Process the pts_df to assign seagrass species based on the seabedtype."""
    pts_df["seagrass_species"] = pts_df["seabedtype"].apply(
        assign_points_seagrass_species
    )
    # Remove undesired cases ("Remove"), drop unused categories
    return pts_df[pts_df["seagrass_species"] != "Remove"]


# Polygon data
def assign_species_habsubtype(habsubtype: str) -> str:
    """Assign seagrass species based on the habsubtype."""
    if habsubtype in ["Cymodocea beds", "Cymodocea nodosa"]:
        return "Cymodocea nodosa"
    elif habsubtype in ["Dead mattes of Posidonia oceanica", "Posidonia oceanica"]:
        return "Posidonia oceanica"
    elif habsubtype in [
        "Zostera beds",
        "Zostera marina",
        "Zostera marina beds",
        "Zostera marina/angustifolia beds",
    ]:
        return "Zostera marina"
    elif habsubtype in ["Zostera noltii", "Zostera noltii beds"]:
        return "Zostera noltei"
    elif habsubtype in ["Unknown"]:
        return "Unspecified"
    elif habsubtype in ["Ruppia beds"]:
        return "Remove"
    else:
        return "Other"


def assign_species_eunis(row: pd.Series) -> str:
    """Assign seagrass species based on the eunis_code."""
    code = row["eunis_code"]
    curr_species = row["seagrass_species"]
    if code in ["A5.531", "A5.5312", "A5.5313", "A5.53131", "A5.53132"]:
        return "Cymodocea nodosa"
    elif code in ["A5.535", "A5.5352", "A5.5353", "A5.5354", "MB252"]:
        return "Posidonia oceanica"
    elif code in ["A5.5331", "A5.533", "A5.5333"]:
        return "Zostera marina"
    elif code in ["A2.61", "A2.611", "A2.6111", "A5.5332"]:
        return "Zostera noltei"
    elif code in ["A5.53", "MB522"]:
        return "Unspecified"
    elif code in ["A5.529", "A5.5343", "A6.61"]:
        return "Remove"
    else:
        return curr_species


def assign_species_anxi(row: pd.Series) -> str:
    """Assign seagrass species based on the anxi_name."""
    anxi_name = row["anxi_name"]
    curr_species = row["seagrass_species"]
    if anxi_name in ["Posidonia beds", "Posidonia oceanica beds"]:
        return "Posidonia oceanica"
    else:
        return curr_species


def process_poly_df_species(poly_df) -> pd.DataFrame:
    """Process the poly_df to assign seagrass species based on the habsubtype, eunis_code, and anxi_name."""
    poly_df["seagrass_species"] = poly_df["habsubtype"].apply(
        assign_species_habsubtype
    )  # correct
    poly_df["seagrass_species"] = poly_df.apply(assign_species_eunis, axis=1)  # correct
    poly_df["seagrass_species"] = poly_df.apply(assign_species_anxi, axis=1)  # correct
    # drop 'remove' values
    return poly_df[poly_df["seagrass_species"] != "Remove"]
