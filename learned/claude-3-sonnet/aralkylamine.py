"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS patterns
    patterns = {
        # Aromatic systems (including heterocycles)
        "aromatic": "[a]1[a][a][a][a]1,[a]1[a][a][a][a][a]1",
        
        # Basic amine patterns
        "primary_amine": "[NX3H2][CX4]",
        "secondary_amine": "[NX3H1;!$(NC=O);!$(N=*);!$(N[a])][CX4]",
        "tertiary_amine": "[NX3H0;!$(NC=O);!$(N=*);!$(N[a])][CX4]",
        
        # Connection patterns
        "aralkyl": "[a]~[CX4]~[CX4;!$(C(=O));!$(C(=N))]~[NX3;!$(NC=O);!$(N=*);!$(N[a])]",
        
        # Exclusion patterns
        "peptide": "[NX3][CX3](=[OX1])[CX4][NX3]",
        "amino_acid": "[NX3H2][CX4H1]([*])[CX3](=[OX1])[OX2H1]",
        "guanidine": "[NX3][CX3](=[NX2])[NX3]"
    }
    
    # Convert patterns to RDKit molecules
    compiled_patterns = {}
    for name, pattern in patterns.items():
        compiled_pattern = Chem.MolFromSmarts(pattern)
        if compiled_pattern is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = compiled_pattern

    # Check for aromatic system
    if not mol.HasSubstructMatch(compiled_patterns["aromatic"]):
        return False, "No aromatic system found"

    # Check for excluded structures
    if mol.HasSubstructMatch(compiled_patterns["peptide"]) or \
       mol.HasSubstructMatch(compiled_patterns["amino_acid"]) or \
       mol.HasSubstructMatch(compiled_patterns["guanidine"]):
        return False, "Contains excluded group (peptide, amino acid, or guanidine)"

    # Check for presence of amine group
    has_amine = False
    for amine_type in ["primary_amine", "secondary_amine", "tertiary_amine"]:
        if mol.HasSubstructMatch(compiled_patterns[amine_type]):
            has_amine = True
            break
    
    if not has_amine:
        return False, "No suitable amine group found"

    # Check for aralkyl connection
    matches = mol.GetSubstructMatches(compiled_patterns["aralkyl"])
    if not matches:
        return False, "No proper connection between aromatic ring and amine"

    # Additional check for molecular properties
    if mol.GetNumAtoms() > 100:  # Avoid very large molecules
        return False, "Molecule too complex"

    return True, "Contains amine group connected to aromatic ring via alkyl linker"