import os
import sys

import PyFBA


def read_media_file(mediaf):
    """
    Read a media file and return a set with the media added. If the environment variable PYFBA_MEDIA_DIR
    is set, we will look in there for mediaf if we can not find it.
        
    Returns a set of compounds that are in the media.

    :param mediaf: The file to read
    :type mediaf: str
    :return: A set of media components
    :rtype: set of metabolism.Compound
    """

    media = set()

    if not os.path.exists(mediaf):
        if 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediaf)):
            mediaf = os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediaf)
        else:
            raise IOError("Media file {} can not be found\nPlease set the environment variable PYFBA_MEDIA_DIR to point to a directory with all the media files".format(mediaf))

    with open(mediaf, 'r') as f:
        for li, l in enumerate(f):
            # skip the header line
            if li == 0:
                continue
            p = l.strip().split("\t")
            if len(p) < 2:
                sys.stderr.write("Skipped line {} as it does not have enough columns\n".format(l.strip()))
                continue
            c = PyFBA.metabolism.Compound(p[1], 'e')
            media.add(c)
    
    return media


def get_pm_plates():
    """
    Build Dictionary of PM plates

    :return: PM plate collections
    :rtype: dict
    """
    plate_media = {"pm_plate_1": [], "pm_plate_2": []}

    # PM PLATE 1
    plate_media["pm_plate_1"] = \
        ["MOPS_NoC_2-Deoxy-D-Ribose.txt",
         "MOPS_NoC_4-Hydroxy-Phenylacetate.txt",
         "MOPS_NoC_Acetic_Acid.txt",
         "MOPS_NoC_Adenosine.txt",
         "MOPS_NoC_Adonitol.txt",
         "MOPS_NoC_Alpha-D-Glucose.txt",
         "MOPS_NoC_Alpha-D-Lactose.txt",
         "MOPS_NoC_Cellobiose.txt",
         "MOPS_NoC_Citric_Acid.txt",
         "MOPS_NoC_D-Alanine.txt",
         "MOPS_NoC_D-Arabinose.txt",
         "MOPS_NoC_D-Arabitol.txt",
         "MOPS_NoC_D-Asparagine.txt",
         "MOPS_NoC_D-Aspartic_Acid.txt",
         "MOPS_NoC_D-Cysteine.txt",
         "MOPS_NoC_D-Fructose.txt",
         "MOPS_NoC_D-Galactose.txt",
         "MOPS_NoC_D-Glucosamine.txt",
         "MOPS_NoC_D-Glucose-6-phosphate.txt",
         "MOPS_NoC_D-Glucose.txt",
         "MOPS_NoC_D-Glutamic_Acid.txt",
         "MOPS_NoC_D-Mannose.txt",
         "MOPS_NoC_D-Ribose.txt",
         "MOPS_NoC_D-Serine.txt",
         "MOPS_NoC_Dulcitol.txt",
         "MOPS_NoC_D-Xylose.txt",
         "MOPS_NoC_Erythritol.txt",
         "MOPS_NoC_Glycerol.txt",
         "MOPS_NoC_Glycine.txt",
         "MOPS_NoC_Inosine.txt",
         "MOPS_NoC_Lactic_Acid.txt",
         "MOPS_NoC_Lactulose.txt",
         "MOPS_NoC_L-Alanine.txt",
         "MOPS_NoC_L-Arabinose.txt",
         "MOPS_NoC_L-Arabitol.txt",
         "MOPS_NoC_L-Asparagine.txt",
         "MOPS_NoC_L-Aspartic_Acid.txt",
         "MOPS_NoC_L-Cysteic_Acid.txt",
         "MOPS_NoC_L-Cysteine.txt",
         "MOPS_NoC_L-Fucose.txt",
         "MOPS_NoC_L-Glutamic_Acid.txt",
         "MOPS_NoC_L-Glutamine.txt",
         "MOPS_NoC_L-Isoleucine.txt",
         "MOPS_NoC_L-Leucine.txt",
         "MOPS_NoC_L-Lysine.txt",
         "MOPS_NoC_L-Methionine.txt",
         "MOPS_NoC_L-Phenylalanine.txt",
         "MOPS_NoC_L-Pyro-Glutamic_Acid.txt",
         "MOPS_NoC_L-Rhamnose.txt",
         "MOPS_NoC_L-Serine.txt",
         "MOPS_NoC_L-Sorbose.txt",
         "MOPS_NoC_L-Threonine.txt",
         "MOPS_NoC_L-Tryptophan.txt",
         "MOPS_NoC_L-Valine.txt",
         "MOPS_NoC_L-Xylose.txt",
         "MOPS_NoC_Malic_Acid.txt",
         "MOPS_NoC_Melibiose.txt",
         "MOPS_NoC_Myo-Inositol.txt",
         "MOPS_NoC_Negative_Control.txt",
         "MOPS_NoC_Oxalic_Acid.txt",
         "MOPS_NoC_Potassium_Sorbate.txt",
         "MOPS_NoC_Propionate.txt",
         "MOPS_NoC_Putrescine.txt",
         "MOPS_NoC_Pyruvate.txt",
         "MOPS_NoC_Quinate.txt",
         "MOPS_NoC_Raffinose.txt",
         "MOPS_NoC_Salicoside.txt",
         "MOPS_NoC_Succinate.txt",
         "MOPS_NoC_Sucrose.txt",
         "MOPS_NoC_Thymidine.txt",
         "MOPS_NoC_Trehalose.txt",
         "MOPS_NoC_Xylitol.txt",
         "MOPS_NoN_Adenine.txt",
         "MOPS_NoN_Adenosine.txt",
         "MOPS_NoN_Allantoin.txt",
         "MOPS_NoN_Beta-Phenylethylamine.txt",
         "MOPS_NoN_Biuret.txt",
         "MOPS_NoN_Cytidine.txt",
         "MOPS_NoN_Cytosine.txt",
         "MOPS_NoN_D-Methionine.txt",
         "MOPS_NoN_D-Valine.txt",
         "MOPS_NoN_Glycine.txt",
         "MOPS_NoN_Guanidine.txt",
         "MOPS_NoN_Histamine.txt",
         "MOPS_NoN_Inosine.txt",
         "MOPS_NoN_L-Arginine.txt",
         "MOPS_NoN_L-Glutathione.txt",
         "MOPS_NoN_L-Histidine.txt",
         "MOPS_NoN_L-Proline.txt",
         "MOPS_NoN_L-Pyro-Glutamic_Acid.txt",
         "MOPS_NoN_N-Acetyl-D-Glucosamine.txt",
         "MOPS_NoN_Negative_Control.txt",
         "MOPS_NoN_Thiourea.txt",
         "MOPS_NoN_Thymine.txt",
         "MOPS_NoN_Tyramine.txt",
         "MOPS_NoN_Uridine.txt"]

    # PM PLATE 2
    plate_media["pm_plate_2"] = \
        ["MOPS_NoN_2-Deoxy-D-Ribose.txt",
         "MOPS_NoN_Acetamide.txt",
         "MOPS_NoN_Ammonium_Chloride.txt",
         "MOPS_NoN_D-Alanine.txt",
         "MOPS_NoN_D-Asparagine.txt",
         "MOPS_NoN_D-Aspartic_Acid.txt",
         "MOPS_NoN_D-Cysteine.txt",
         "MOPS_NoN_D-Glucosamine.txt",
         "MOPS_NoN_D-Glutamic_Acid.txt",
         "MOPS_NoN_D-Serine.txt",
         "MOPS_NoN_L-Alanine.txt",
         "MOPS_NoN_L-Asparagine.txt",
         "MOPS_NoN_L-Citrulline.txt",
         "MOPS_NoN_L-Cysteine.txt",
         "MOPS_NoN_L-Glutamic_Acid.txt",
         "MOPS_NoN_L-Glutamine.txt",
         "MOPS_NoN_L-Isoleucine.txt",
         "MOPS_NoN_L-Leucine.txt",
         "MOPS_NoN_L-Lysine.txt",
         "MOPS_NoN_L-Methionine.txt",
         "MOPS_NoN_L-Ornithine.txt",
         "MOPS_NoN_L-Phenylalanine.txt",
         "MOPS_NoN_L-Serine.txt",
         "MOPS_NoN_L-Threonine.txt",
         "MOPS_NoN_L-Tryptophan.txt",
         "MOPS_NoN_L-Valine.txt",
         "MOPS_NoN_Negative_Control.txt",
         "MOPS_NoN_Putrescine.txt",
         "MOPS_NoN_Tyrosine.txt",
         "MOPS_NoP_Adenosine-5-Monophosphate.txt",
         "MOPS_NoP_Beta-Glycerophosphate.txt",
         "MOPS_NoP_Creatinephosphate.txt",
         "MOPS_NoP_D-Glucose-6-Phosphate.txt",
         "MOPS_NoP_DL-Alpha-Glycerophosphate.txt",
         "MOPS_NoP_Diethyl-Dithiophosphate.txt",
         "MOPS_NoP_Negative_Control.txt",
         "MOPS_NoP_Potassium_Phosphate.txt",
         "MOPS_NoP_Sodium_Pyrophosphate_INCOMPLETE.txt",
         "MOPS_NoP_Sodium_Thiophosphate.txt",
         "MOPS_NoS_1-Butane-Sulfonic_Acid.txt",
         "MOPS_NoS_Acetyl_Cysteine.txt",
         "MOPS_NoS_D-Cysteine.txt",
         "MOPS_NoS_D-Methionine.txt",
         "MOPS_NoS_DL-Alpha-Amino-N-Butyric_Acid_INCOMPLETE.txt",
         "MOPS_NoS_DL-Ethionine.txt",
         "MOPS_NoS_Diethyl-Dithiophosphate_INCOMPLETE.txt",
         "MOPS_NoS_Gamma-Amino-N-Butyric_Acid.txt",
         "MOPS_NoS_Glutathione.txt",
         "MOPS_NoS_Isethionic_Acid.txt",
         "MOPS_NoS_L-Cysteic_Acid.txt",
         "MOPS_NoS_L-Cysteine.txt",
         "MOPS_NoS_L-Djenkolic_Acid.txt",
         "MOPS_NoS_L-Methionine.txt",
         "MOPS_NoS_Magnesium_Sulfate.txt",
         "MOPS_NoS_Methane_Sulfonic_Acid.txt",
         "MOPS_NoS_N-Acetyl-DL-Methionine.txt",
         "MOPS_NoS_N-Acetyl-L-Cysteine.txt",
         "MOPS_NoS_Negative_Control.txt",
         "MOPS_NoS_Potassium-Tetra-Thionate.txt",
         "MOPS_NoS_Sodium_Thiosulfate.txt",
         "MOPS_NoS_Sulfanic_Acid.txt",
         "MOPS_NoS_Taurine.txt",
         "MOPS_NoS_Taurocholic_Acid.txt",
         "MOPS_NoS_Thiourea.txt"]

    return plate_media
