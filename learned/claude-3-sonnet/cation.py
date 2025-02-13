"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:49637 cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for atoms with positive formal charges
    cationic_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0]
    if cationic_atoms:
        # Exclude nitro groups and metal complexes
        if not any(atom.GetIsAromatic() and atom.GetAtomicNum() == 7 and sum(nbo.GetBondOrder() for nbo in atom.GetBonds()) == 3 for atom in mol.GetAtoms()):
            if not any(atom.GetAtomicNum() > 20 for atom in mol.GetAtoms()):
                return True, "Contains atoms with positive formal charges"
    
    # Check for bracketed cationic elements
    cationic_elements = Chem.MolFromSmarts("[+]")
    cationic_element_matches = mol.GetSubstructMatches(cationic_elements)
    if cationic_element_matches:
        return True, "Contains bracketed cationic elements"
    
    # Check for common organic cation substructures
    ammonium_pattern = Chem.MolFromSmarts("[NH3+]")
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+]")
    sulfonium_pattern = Chem.MolFromSmarts("[S+]")
    phosphonium_pattern = Chem.MolFromSmarts("[P+]")
    oxonium_pattern = Chem.MolFromSmarts("[O+]")
    
    organic_cation_patterns = [ammonium_pattern, quaternary_ammonium_pattern, sulfonium_pattern, phosphonium_pattern, oxonium_pattern]
    for pattern in organic_cation_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {pattern.GetSmarts()} substructure"
    
    # Check for zwitterionic species
    zwitterion_query = rdqueries.HasChargeDistributionMatchingDescription("Zwitterion")
    if zwitterion_query.IsSatisfied(mol):
        return True, "Molecule is a zwitterion"
    
    return False, "No cationic features found"