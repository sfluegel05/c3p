"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    This function uses a combination of substructure searches and property checks.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function to count carbons and oxygens
    def count_atoms(molecule, atom_symbol):
        return sum(1 for atom in molecule.GetAtoms() if atom.GetSymbol() == atom_symbol)

    # Minimum Carbon and oxygen atoms check (to remove small/unlikely molecules)
    min_carbons = 10
    min_oxygens = 1
    c_count = count_atoms(mol, "C")
    o_count = count_atoms(mol, "O")
    if c_count < min_carbons or o_count < min_oxygens:
      return False, f"Too few carbon ({c_count}) or oxygen atoms ({o_count}). Must have at least {min_carbons} C and {min_oxygens} O atoms."


    # 1. Check for a fused ring system containing a 6-membered ring (core of cannabinoids)
    fused_ring_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]2[c]1[c][c][c][c]2")
    if mol.HasSubstructMatch(fused_ring_pattern):
        return True, "Contains a fused ring system characteristic of many cannabinoids"


    # 2. Check for long carbon chains with amides/esters/ethers
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if mol.HasSubstructMatch(long_chain_pattern):
      #check if it has any of the functional groups of interest
        amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2]")
        ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
        ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
        
        if mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(ether_pattern):
            return True, "Contains a long carbon chain with a functional group (amide/ester/ether)."

    #3. Check for indoles or other nitrogen heterocycles.
    indole_pattern = Chem.MolFromSmarts("c1ccnc2c1cc[cH]2")
    if mol.HasSubstructMatch(indole_pattern):
         return True, "Contains an indole heterocycle"
    
    other_n_het_pattern = Chem.MolFromSmarts("N1[c][c][c][c][c]1")
    if mol.HasSubstructMatch(other_n_het_pattern):
        return True, "Contains a nitrogen heterocycle."

    return False, "Does not fit the criteria for a cannabinoid structure."