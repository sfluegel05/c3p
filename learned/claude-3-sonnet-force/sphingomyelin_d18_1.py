"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: CHEBI:17706 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_sphingomyelin_d18_1(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 is a sphingomyelin containing a sphingosine backbone and a C18:1 fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingosine backbone features
    # 1. Check for long aliphatic chain (16-24 carbon atoms)
    longest_chain = Chem.GetLongestAliphaticChain(mol)
    if not (16 <= len(longest_chain) <= 24):
        return False, "Sphingosine chain length not in expected range (16-24 carbons)"

    # 2. Check for presence of a double bond in the chain
    if not any(bond.GetIsAromatic() for bond in longest_chain.GetBonds()):
        return False, "No double bond found in sphingosine chain"

    # 3. Check for presence of a primary hydroxy group (-OH) in the chain
    has_primary_alcohol = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetTotalNumHs() == 1:
            has_primary_alcohol = True
            break
    if not has_primary_alcohol:
        return False, "No primary alcohol group found in sphingosine chain"

    # Look for C18:1 fatty acyl chain
    c18_1_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCC(=O)")
    c18_1_matches = mol.GetSubstructMatches(c18_1_pattern)
    if len(c18_1_matches) != 1:
        return False, "No C18:1 fatty acyl chain found"

    # Check for amide bond between sphingosine and fatty acyl chain
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, "Missing amide bond between sphingosine and fatty acyl chain"

    # Check for phosphocholine head group
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine head group found"

    return True, "Contains sphingosine backbone and C18:1 fatty acyl chain linked via amide bond"