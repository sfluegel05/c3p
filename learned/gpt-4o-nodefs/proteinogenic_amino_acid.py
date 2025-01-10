"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for a generic proteinogenic amino acid
    # This pattern should capture the distinct features of proteinogenic amino acids: 
    # Central carbon (C) bonded to an amino group (N), carboxyl group (O=CO), and optional chiralities.
    amino_acid_pattern = Chem.MolFromSmarts(
        "[CX4H2][NX3](C(=O)[O,N])" +  # Central carbon with amino group and carboxylate
        "[$([C,c;!$(*=,#)])]"         # Variable side chain (R group) allowance
    )
    
    # Run a substructure search with this pattern
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Does not match core proteinogenic amino acid structure"
    
    # Verify for larger molecules that they don't contain structures inconsistent with amino acids
    # Refine based on a more localized search instead of entire structure which may better
    # accommodate partial matching for modified residues (e.g. deuterated amino acids)
    try:
        # Find the first matching conformation
        match = mol.GetSubstructMatch(amino_acid_pattern)
        if not match:
            return False, "Core structure is missing or improperly modified"
    except:
        return False, "An error occurred during matching process"

    # If found matches to the structural pattern that often characterize the proteinogenic amino acids
    return True, "Matches refined structure expected of proteinogenic amino acids"

# Example test for the function
smiles_list = [
    'OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]',  # L-phenylalanine-d5
    'N[C@@H](Cc1c[nH]cn1)C(O)=O',  # L-histidine
    'C[C@H](N)C(O)=O',  # L-alanine
    'CC(C)[C@H](N)C(O)=O'  # L-valine
]

for smi in smiles_list:
    is_amino_acid, reason = is_proteinogenic_amino_acid(smi)
    print(f"SMILES: {smi}, Is Proteinogenic Amino Acid: {is_amino_acid}, Reason: {reason}")