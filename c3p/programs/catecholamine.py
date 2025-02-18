"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts
from rdkit.Chem.rdmolops import GetLargestFragment

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine has a benzene ring with two adjacent hydroxyl groups (catechol)
    and a 2-aminoethyl group (-CH2-CH2-NR2) or its substituted derivatives attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split into components and take the largest fragment (to ignore salts, counterions)
    mol = GetLargestFragment(mol)
    
    # Check for benzene ring with two adjacent hydroxyl groups (catechol)
    # SMARTS: two adjacent aromatic hydroxyls on a benzene ring
    catechol_pattern = MolFromSmarts("[OH]c1:c:c:c(:c:c:1)-[OH]")
    if not mol.HasSubstructMatch(catechol_pattern):
        # Try alternative positions for the hydroxyls
        catechol_pattern2 = MolFromSmarts("c1c(O)c(O)ccc1")
        if not mol.HasSubstructMatch(catechol_pattern2):
            return False, "No catechol group (adjacent aromatic hydroxyls) found"
    
    # Check for 2-aminoethyl group or substituted derivatives directly attached to the ring
    # SMARTS: aromatic carbon connected to two carbons (CH2) followed by an amine
    ethylamine_pattern = MolFromSmarts("[c;r6]-[CH2]-[CH2]-[NX3;H0,H1,H2;!$(N=*)]")
    if not mol.HasSubstructMatch(ethylamine_pattern):
        # Allow branching on the ethyl chain (e.g., -CH(CH3)-NH2)
        ethylamine_pattern2 = MolFromSmarts("[c;r6]-[CH2]-[CX4]-[NX3;H0,H1,H2;!$(N=*)]")
        if not mol.HasSubstructMatch(ethylamine_pattern2):
            # Allow substituents on the amine (e.g., -NHCH3)
            ethylamine_pattern3 = MolFromSmarts("[c;r6]-[CX4]-[CX4]-[NX3;H0,H1,H2;!$(N=*)]")
            if not mol.HasSubstructMatch(ethylamine_pattern3):
                return False, "No 2-aminoethyl-like group attached to aromatic ring"
    
    return True, "Contains catechol group with 2-aminoethyl-like substituent"