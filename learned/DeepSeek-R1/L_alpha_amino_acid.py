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
    Must have an amino group and carboxyl group on adjacent carbons (alpha-carbon),
    with L-configuration (S CIP code) at the alpha-carbon.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all carboxyl groups (COOH or COO-)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"
    
    # Check each carboxyl group's adjacent carbon (alpha-carbon)
    for carboxyl in carboxyl_matches:
        carb_carbon = carboxyl[0]  # The carbon in the carboxyl group
        # Get the alpha-carbon (connected to carboxyl carbon)
        alpha_carbons = [n.GetIdx() for n in mol.GetAtomWithIdx(carb_carbon).GetNeighbors() if n.GetAtomicNum() == 6]
        for alpha_idx in alpha_carbons:
            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            # Check if alpha-carbon has an amino group (NH2, NH, etc.)
            amino_ns = [n for n in alpha_atom.GetNeighbors() if n.GetAtomicNum() == 7]
            if not amino_ns:
                continue
            # Check if any of the nitrogens is an amine (at least one H or part of a ring)
            has_amino = False
            for n in amino_ns:
                if n.GetTotalNumHs() > 0 or n.IsInRing():
                    has_amino = True
                    break
            if not has_amino:
                continue
            
            # Check chirality of alpha-carbon
            if alpha_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            
            # Assign CIP configuration
            try:
                Chem.AssignCIPLabels(mol)
            except:
                continue  # Skip if CIP assignment fails
            
            if not alpha_atom.HasProp("_CIPCode"):
                continue
            cip = alpha_atom.GetProp("_CIPCode")
            if cip.upper() == "S":
                return True, "L-alpha-amino acid with S configuration at alpha-carbon"
    
    return False, "No L-alpha-amino acid structure found"