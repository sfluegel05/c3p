"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a specific core structure featuring a furan ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Attempt to parse the SMILES string into a molecular structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sufficient oxygenation, considering known examples have high oxygen counts
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 10:  # Limonoids show significant oxygen content, refine as needed
        return False, "Insufficient oxygen content typical of limonoids"
    
    # Define SMARTS pattern for the furan ring as a common feature
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if furan_pattern is not None and not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring identified, atypical for limonoid"
    
    # Hypothetical SMARTS pattern for a limonoid core structure
    # Needs refinement based on structural diversity of actual limonoids
    limonoid_core_pattern = Chem.MolFromSmarts("C1=CCC2(O)C3=C(O)C(C)C4C(=O)CCC5=CC6C7C1C3=C4C(=O)C(O)C5(C7C6=O)=O")
    if limonoid_core_pattern is not None and not mol.HasSubstructMatch(limonoid_core_pattern):
        return False, "Lacks prototypical limonoid steroid-like core structure"

    # Consider molecular weight based on highly oxygenated limonoid nature
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 450 or mol_wt > 1200:
        return False, "Molecular weight falls outside typical limonoid range"

    return True, "Molecule exhibits characteristics consistent with limonoid structural class"