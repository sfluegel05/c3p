"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: An acyl-CoA is a thioester that results from the condensation of the thiol group 
of Coenzyme A with the carboxy group of a fatty acid. Its structure features an acyl thioester 
group directly attached to a pantetheine moiety, which in turn is connected to an adenosine 
nucleotide.
This improved classifier demands that the thioester is directly linked to the CoA “core”:
   – An acyl carbonyl (C(=O)) attached to a sulfur.
   – The sulfur immediately followed by the pantetheine fragment "SCCNC(=O)CCNC(=O)".
   – The presence of an adenine moiety "n1cnc2c(N)ncnc12" (found in virtually all acyl-CoAs).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    The function requires that the molecule contains:
      1. A thioester group directly connected with the pantetheine linkage. 
         This is captured by the SMARTS "[CX3](=O)[SX2]CCNC(=O)CCNC(=O)".
      2. The adenosine (adenine) moiety attached to the ribose‐phosphate in CoA,
         detected by the SMARTS "n1cnc2c(N)ncnc12".
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define a SMARTS pattern for the acyl thioester directly linked to the pantetheine part.
    # This pattern means: a carbonyl C with a double-bond O attached to a sulfur that is immediately
    # followed by the sequence: CCNC(=O)CCNC(=O)
    acyl_coa_core_smarts = "[CX3](=O)[SX2]CCNC(=O)CCNC(=O)"
    acyl_coa_core = Chem.MolFromSmarts(acyl_coa_core_smarts)
    if acyl_coa_core is None:
        return False, "Error creating acyl-CoA core SMARTS pattern"
    
    if not mol.HasSubstructMatch(acyl_coa_core):
        return False, "Acyl-CoA core not found (required thioester linked to pantetheine missing)"

    # 2. Define a SMARTS pattern for the adenine moiety of Coenzyme A.
    adenine_smarts = "n1cnc2c(N)ncnc12"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if adenine is None:
        return False, "Error creating adenine SMARTS pattern"
    
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety of CoA not detected"
    
    # If both key fragments are identified, we consider the molecule as an acyl-CoA.
    return True, "Contains acyl thioester linked to CoA pantetheine and adenine moiety"

# Example usage:
# status, reason = is_acyl_CoA("[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)...")
# print(status, reason)