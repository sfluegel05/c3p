"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives
of these compounds.
Improved classifier that relaxes the molecular weight threshold and requires a long alkyl chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds include sphinganine, its homologs/stereoisomers, and its hydroxy/unsaturated derivatives.
    They typically possess a long aliphatic chain and a characteristic amino-diol (or related carbonyl) core.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as sphingoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a long aliphatic chain.
    # We require a contiguous chain of at least 8 sp3 carbons.
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long aliphatic chain (>=8 carbons) found."

    # Count carbon atoms (a rough measure -- many sphingoid molecules have many C atoms).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Too few carbon atoms ({c_count}); sphingoid molecules generally contain a long chain."

    # Check that at least one nitrogen is present, as sphingoid bases have an amino group.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found; sphingoid compounds require an amine group in their backbone."
    
    # Define SMARTS patterns for the core sphingoid motif.
    # Pattern 1: typical amino-diol: C(O)C(N)CO (ignoring stereochemistry)
    core_pattern1 = Chem.MolFromSmarts("C(O)C(N)CO")
    # Pattern 2: a carbonyl variant (observed in unsaturated or dehydro forms): C(=O)C(N)CO
    core_pattern2 = Chem.MolFromSmarts("C(=O)C(N)CO")
    
    core_found = mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)
    if not core_found:
        return False, "No typical sphingoid core pattern (amino-diol or carbonyl variant) found."
    
    # Check molecular weight: many sphingoid bases are around 280-350 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical sphingoid compound."
    
    return True, "Molecule exhibits a sphingoid core with a long aliphatic chain and appropriate functionality."

# Example test cases:
if __name__ == "__main__":
    example_smiles = [
        "CCCCCCCCCCCC\\C=C\\C(=O)[C@@H](N)CO",  # 3-dehydrosphingosine (should return True)
        "OC[C@@]([C@@](CCCCCCCCCCCCCC)(O)H)(N)H",  # heptadecasphinganine (should return True)
        "CCO",  # too short, not sphingoid
    ]
    for smi in example_smiles:
        res, reason = is_sphingoid(smi)
        print(f"SMILES: {smi}\n Classified: {res}\n Reason: {reason}\n")