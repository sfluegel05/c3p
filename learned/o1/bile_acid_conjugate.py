"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a hydrophilic group such as glycine,
    taurine, amino acids, sulfate, glucuronic acid, glucose, other sugars, or coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a general pattern for the steroid nucleus (four fused rings)
    steroid_pattern = Chem.MolFromSmarts("""
    [#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4
    """)
    
    # Define a pattern for the bile acid side chain ending with carboxylic acid at C24
    side_chain_pattern = Chem.MolFromSmarts("""
    [#6]-[#6]-[#6]-[#6]-C(=O)[O;H,–]
    """)
    
    # Combine steroid nucleus and side chain to define bile acid core
    bile_acid_pattern = Chem.MolFromSmarts("""
    [$([#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4),R]
    -[*]-[*]-[*]-C(=O)[O;H,–]
    """)
    
    # Check for bile acid core
    if not mol.HasSubstructMatch(bile_acid_pattern):
        return False, "No bile acid core found"

    # Define SMARTS patterns for conjugated groups attached via amide or ester linkage
    conjugates = {
        "glycine": Chem.MolFromSmarts("C(=O)NCC(=O)[O;H]"),
        "taurine": Chem.MolFromSmarts("C(=O)NCCS(=O)(=O)[O;H]"),
        "amino_acid": Chem.MolFromSmarts("C(=O)N[C@@H]([#6])[#6]C(=O)[O;H]"),  # General amino acid
        "sulfate": Chem.MolFromSmarts("OS(=O)(=O)[O;H]"),
        "glucuronic_acid": Chem.MolFromSmarts("""
        [#6]-1(-[#8]-[#6]-[#8]C(=O)[O;H])[#6]([#8])-[#6]([#8])-[#6]([#8])-[#6]-1
        """),
        "glucose": Chem.MolFromSmarts("""
        [#6]-1(-[#8]-[#6])[#6]([#8])-[#6]([#8])-[#6]([#8])-[#6]-1
        """),
        "coenzyme_A": Chem.MolFromSmarts("NC(=O)CCNC(=O)CC(=O)NCCS"),
    }
    
    # Check for conjugation with any of the specified groups
    conjugated = False
    conjugate_name = ""
    for name, pattern in conjugates.items():
        # Attachment via an amide or ester bond to the bile acid carboxyl group
        linker_pattern = Chem.MolFromSmarts(f"C(=O)[O,N]*{pattern}")
        if mol.HasSubstructMatch(linker_pattern):
            conjugated = True
            conjugate_name = name
            break

    if not conjugated:
        return False, "No conjugated hydrophilic group found"

    return True, f"Contains bile acid core conjugated with {conjugate_name}"