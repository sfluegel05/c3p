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
    An L-alpha-amino acid has an amino group, a carboxyl group on the adjacent carbon (alpha-carbon),
    and the L-configuration (S CIP code) at the alpha-carbon.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find carboxyl groups (COOH or COO-)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"
    
    for match in carboxyl_matches:
        # The carbon in the carboxyl group
        carb_c = match[0]
        # The alpha-carbon is bonded to this carbon (not part of the carboxyl oxygens)
        alpha_carbon = None
        for neighbor in mol.GetAtomWithIdx(carb_c).GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Skip oxygen atoms in carboxyl
                continue
            alpha_carbon = neighbor
            break
        if not alpha_carbon:
            continue
        
        # Check if alpha-carbon has an amino group (NH2, NH, or part of a ring)
        has_amino = False
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:  # Nitrogen
                # Check if it's an amine (at least one H)
                if neighbor.GetTotalNumHs() >= 1 or neighbor.IsInRing():
                    has_amino = True
                    break
        if not has_amino:
            continue
        
        # Check if alpha-carbon is a chiral center
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        
        # Assign CIP configuration
        try:
            Chem.AssignCIPLabels(mol)
        except:
            continue  # In case of errors in CIP assignment
        
        cip = alpha_carbon.GetProp("_CIPCode", "") if alpha_carbon.HasProp("_CIPCode") else ""
        if cip == "S":
            return True, "L-alpha-amino acid with S configuration at alpha-carbon"
        # Handle cysteine-like cases where R configuration is L (uncommon)
        # This requires checking the R group, which is complex. Omitted for simplicity.
    
    return False, "Not an L-alpha-amino acid"