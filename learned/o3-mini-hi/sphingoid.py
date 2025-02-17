"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives
of these compounds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds include sphinganine, its homologs/stereoisomers, and its hydroxy/unsaturated derivatives.
    They typically have a long aliphatic chain and a characteristic amino-diol (or related carbonyl) core.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as sphingoid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbon atoms ({c_count}); sphingoid molecules generally contain a long aliphatic chain."
    
    # Ensure the molecule has at least one nitrogen (sphingoid bases have a primary/secondary amine)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found; sphingoid compounds require an amine group in their backbone."
    
    # Define SMARTS patterns for the sphingoid core.
    # Pattern 1: a typical amino-diol motif: C(O)C(N)CO 
    core_pattern1 = Chem.MolFromSmarts("C(O)C(N)CO")
    # Pattern 2: a carbonyl variant found in dehydro forms: C(=O)C(N)CO
    core_pattern2 = Chem.MolFromSmarts("C(=O)C(N)CO")
    
    core_match1 = mol.HasSubstructMatch(core_pattern1)
    core_match2 = mol.HasSubstructMatch(core_pattern2)
    if not (core_match1 or core_match2):
        return False, "No typical sphingoid core pattern (amino-diol or related motif) found."
    
    # (Optional additional check: molecular weight should be within a range
    # expected for sphingoid bases; many examples are above ~300 Da.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical sphingoid compound."
    
    return True, "Molecule exhibits a sphingoid core with a long aliphatic chain and appropriate functionality."

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCC\C=C\[C@@H](O)[C@@H](N)CO"  # 3-dehydrosphingosine example
# result, reason = is_sphingoid(test_smiles)
# print(result, reason)