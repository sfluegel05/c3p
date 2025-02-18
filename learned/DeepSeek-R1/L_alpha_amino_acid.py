"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:32644 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    Must have an amino group and carboxyl group on the alpha-carbon with S configuration (L-form).
    Excludes peptide-bonded amino acids.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for alpha-amino acid core: NH2-CH(R)-COOH
    # Allows for variations in protonation and substitution
    amino_acid_core = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]-[C]([CX3](=O)[OX2H1,O-])-*")
    core_matches = mol.GetSubstructMatches(amino_acid_core)
    
    if not core_matches:
        return False, "No alpha-amino acid core found"
    
    # Check each potential alpha-carbon for L-configuration
    for match in core_matches:
        alpha_carbon = match[1]  # Index 1 is the central carbon (alpha)
        alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
        
        # Verify chirality exists
        if alpha_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        
        # Assign CIP configuration
        try:
            Chem.AssignCIPLabels(mol)
        except:
            continue
        
        if not alpha_atom.HasProp("_CIPCode"):
            continue
        
        cip = alpha_atom.GetProp("_CIPCode")
        if cip.upper() == "S":
            return True, f"L-alpha-amino acid with S configuration at alpha-carbon (atom {alpha_carbon})"
    
    return False, "No L-alpha-amino acid structure with S configuration found"