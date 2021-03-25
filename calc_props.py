#
# Calculation of molecular descriptors and complexity index
#
# Implements RDKit
#
# Bryon Drown, May 2015
# Updated Oct. 9, 2015
# University of Illinois, Urbana-Champaign
#
__doc__ = """
Performs calculations of physiochemical properties of set of compounds
Properties to be calculated:
Fsp3, chiral centers
"""

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from collections import defaultdict
from collections import OrderedDict
import optparse
from os.path import basename
import csv


def main():

    args = parse_args()

    ms = [x for x in Chem.SDMolSupplier(args.input_file) if x is not None]

    mols = []
    for mol in ms:
        remove_ligprep_props(mol)
        temp = calc_builtin_props(mol)
        calcRingDescriptors(temp)
        calc_chiral_centers(temp)
        mols.append(temp)

    filename_base = args.output_path + '/' + \
        basename(args.input_file).strip((".sdf"))
    if(args.csv):
        with open(filename_base + '.csv', 'w') as out:
            write_mol_csv(mols, out)
    if(args.sdf):
        write_mol_sdf(mols, filename_base)


def write_mol_csv(mols, outfile, includeChirality=True):
    """Writes list of molecules and properties to CSV file
    """
    w = csv.writer(outfile)

    # Get prop names from first mol and get header
    first = mols[0]
    propNames = list(first.GetPropNames())
    outL = []
    outL.append('SMILES')
    outL.extend(propNames)
    w.writerow(outL)

    # Write out properties for each molecule
    for mol in mols:
        smi = Chem.MolToSmiles(mol, isomericSmiles=includeChirality)
        outL = []
        outL.append(smi)
        for prop in propNames:
            if mol.HasProp(prop):
                outL.append(str(mol.GetProp(prop)))
            else:
                outL.append('')
        w.writerow(outL)
    return


def csv_header(mol):
    """Creates string for header of csv file
    """
    properties = ['name', 'smiles']
    props = mol.GetPropNames()
    for prop in props:
        properties.append(prop)
    return ','.join(map(str, properties))


def mol_props_to_csv(mol):
    """Creates string for properties of individual molecule when writing to csv
    TODO: keep width the same when some compounds don't have a given property
    """
    values = []
    values.append(mol.GetProp('_Name'))
    values.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    properties = mol.GetPropNames()
    for prop in properties:
        values.append(mol.GetProp(prop))

    return ','.join(map(str, values))


def write_mol_sdf(mols, filename):
    """Writes list of molecules and properties to SDF file
    """
    ms_wr = Chem.SDWriter(filename + ".sdf")
    for mol in mols:
        ms_wr.write(mol)
    ms_wr.close()


def parse_args():
    """Parse the command line options.

    @return:  All script options
    """

    parser = optparse.OptionParser(__doc__)
    parser.set_defaults(verbose=False)
    parser.add_option("-i", "--input", dest="input_file", default=None,
                      help="Input sdf file that contains structures for which properties will be calculated [default: %default]")
    parser.add_option("-o", "--output", dest="output_path", default='',
                      help="Folder to which output files should be saved [default: %default]")
    parser.add_option("-c", "--csv", action="store_true", dest="csv")
    parser.add_option("-s", "--sdf", action="store_true", dest="sdf")

    (options, args) = parser.parse_args()
    if(options.input_file == None):
        print("Input file %s is needed" % options.input_file)
        parser.print_help()
        sys.exit(1)
    return options


def remove_ligprep_props(m):
    """The properties that are attached to molecules by ligprep were are removed
    """
    props = ['chiral flag', 'version', 's_m_source_file', 'i_m_source_file_index', 'i_lp_mmshare_version', 'r_lp_tautomer_probability',
             'r_epik_Ionization_Penalty', 'r_epik_Ionization_Penalty_Charging', 'r_epik_Ionization_Penalty_Neutral',
             'r_epik_State_Penalty', 'r_epik_Charging_Adjusted_Penalty', 'i_epik_Tot_Q', 'i_epik_Tot_abs_Q', 'i_f3d_flags',
             's_lp_Force_Field', 'r_lp_Energy', 'b_lp_Chiralities_Consistent', 's_lp_Variant', 's_epik_Chemistry_Notes', 's_epik_input',
             's_epik_cmdline']
    for prop in props:
        m.ClearProp(prop)


def calc_builtin_props(m):
    """Calculates properties that are part of rdkit base

    @param m: molecule for which to perform calculations
    @return: molecule with properties attached
    """

    nms = ('FractionCSP3', 'MolWt', 'RingCount')
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)

    descrs = calc.CalcDescriptors(m)
    for x in range(len(descrs)):
        m.SetProp(str(nms[x]), str(descrs[x]))

    return m

def calc_chiral_centers(m):
    """Calculates the number of chiral centers in a molecule

    @param m: molecule for which to perform calculations
    """

    centers = Chem.FindMolChiralCenters(m, force=True, includeUnassigned=True)
    m.SetProp('NumChiralCenters', str(len(centers)))
    return



def calcRingDescriptors(m):
    """Calculates a set of properties that measure ring complexity

    @param m: molecule for which to perform calculations
    """

    nBonds = m.GetNumBonds()
    nAtoms = m.GetNumAtoms()
    cyclomatic = nBonds - nAtoms + 1
    if(cyclomatic < 1):
        return

    ri = m.GetRingInfo()
    if(ri.NumRings() < 1):
        return
    # get total ring path and nBondRings
    totalRing = 0
    Bonds = []
    Bridges = []
    for ring in ri.BondRings():

        for id in ring:

            if (ri.NumBondRings(id) > 1):

                Bridges.append(id)
            totalRing += 1
            Bonds.append(id)

    # remove duplicates, then get length
    nBondRings = len(OrderedDict.fromkeys(Bonds).keys())
    nBridgeEdges = len(OrderedDict.fromkeys(Bridges).keys())

    # get nAtomRings
    Atoms = []
    for ring in ri.AtomRings():

        for id in ring:

            Atoms.append(id)
    nAtomRings = len(OrderedDict.fromkeys(Atoms).keys())

    # descriptors
    ringFusionDensity = 2 * float(nBridgeEdges) / float(nAtomRings)
    ringComplexityIndex = float(totalRing) / float(nAtomRings)
    molecularCyclizedDegree = float(nAtomRings) / float(nAtoms)
    nRingSystems = (nBonds - nBondRings) - (nAtoms - nAtomRings) + 1
    if(nRingSystems < 1):

        ringFusionDegree = 0
    else:

        ringFusionDegree = float(cyclomatic) / float(nRingSystems)

    # set props
    m.SetProp('TotalRing', str(totalRing))
    m.SetProp('NumBridges', str(nBridgeEdges))
    m.SetProp('nBondRings', str(nBondRings))
    m.SetProp('nAtomRings', str(nAtomRings))
    m.SetProp('ringFusionDensity', str(ringFusionDensity))
    m.SetProp('ringFusionDegree', str(ringFusionDegree))
    m.SetProp('ringComplexityIndex', str(ringComplexityIndex))
    m.SetProp('molecularCyclizedDegree', str(molecularCyclizedDegree))
    m.SetProp('NumRingSystems', str(nRingSystems))

    return


if __name__ == '__main__':
    main()
