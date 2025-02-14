"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Complexity: Number of rings and chiral centers
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    if num_rings < 1 and num_chiral_centers < 1:
         return False, "Too simple, lacks ring structures and chiral centers"

    if num_heavy_atoms < 15:
        return False, "Too few heavy atoms to be a semisynthetic derivative."

    # 2. Look for common synthetic modifications

    # Ester pattern
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    # Amide pattern
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    num_amides = len(mol.GetSubstructMatches(amide_pattern))
    # Ether pattern
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    # Halogens
    num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])


    if num_esters + num_amides + num_ethers == 0 and num_halogens == 0 :
            return False, "No common synthetic modifications found"


    #3. Exclude simple molecules (already done above)

    #4. Exclude fully synthetic compounds
    #Check for many halogens and lack of chiral centers
    if num_halogens > 2 and num_chiral_centers < 2:
        return False, "Likely a fully synthetic molecule with many halogens and few chiral centers."


    return True, "Contains structural features suggestive of a semisynthetic derivative"