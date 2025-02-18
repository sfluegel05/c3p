"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
Definition: The amino acid ester derivative obtained by the formal condensation of 
an alpha-amino acid with an alcohol.

An alpha-amino acid typically has the structure NH2–CHR–COOH. Upon esterification 
the acid (COOH) is converted to an ester (COOR). Here we look for a substructure
where an amino group is connected to an sp3 carbon (the α‐carbon) which is directly 
bound to a carboxyl group converted into an ester (C(=O)O[*]).
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    It searches for a pathway corresponding to the motif: [NX3]-[CX4]-[C](=O)O[*],
    which corresponds to an amino acid (NH2–CHR–COOH) that has been esterified.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an alpha-amino acid ester motif, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern corresponding to N–alphaC–C(=O)O[*]
    # This pattern looks for an amino group [NX3] attached to an sp3 carbon ([CX4]),
    # which is connected to a carbonyl carbon [C](=O) bonded to an oxygen that is bound to any atom ([*]).
    aa_ester_smarts = "[NX3]-[CX4]-[C](=O)O[*]"
    aa_ester_pattern = Chem.MolFromSmarts(aa_ester_smarts)
    if aa_ester_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Check if the molecule has at least one match to the alpha-amino acid ester motif
    if not mol.HasSubstructMatch(aa_ester_pattern):
        return False, "No alpha-amino acid ester motif found"
    
    return True, "Contains an alpha-amino acid ester moiety"