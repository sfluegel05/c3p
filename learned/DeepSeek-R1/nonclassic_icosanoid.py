"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are biologically active C20 fatty acid derivatives with oxygenation, excluding prostanoids and leukotrienes.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for carboxylic acid or ester (essential for fatty acid derivatives)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    has_acid = mol.HasSubstructMatch(acid_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    if not (has_acid or has_ester):
        return False, "No carboxylic acid/ester group"

    # Collect positions of acid/ester carbonyl carbons for exclusion
    acid_ester_carbonyls = set()
    for match in mol.GetSubstructMatches(acid_pattern) + mol.GetSubstructMatches(ester_pattern):
        acid_ester_carbonyls.add(match[0])

    # Oxygenation checks (must have at least one beyond acid/ester)
    has_oxygenation = False
    # Check hydroxyls not in acid group
    acid_oxygens = {match[1] for match in mol.GetSubstructMatches(acid_pattern)}
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) and any(oh[0] not in acid_oxygens for oh in mol.GetSubstructMatches("[OH]")):
        has_oxygenation = True
    # Check epoxy groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")):
        has_oxygenation = True
    # Check ketones not in acid/ester
    if any(match[0] not in acid_ester_carbonyls for match in mol.GetSubstructMatches("[CX3]=[OX1]")):
        has_oxygenation = True

    if not has_oxygenation:
        return False, "No hydroxyl/epoxy/ketone groups beyond acid/ester"

    # Exclude prostanoids (5-membered all-carbon rings)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            return False, "Contains prostanoid ring"

    # Check approximate C20 backbone via chain length from carboxylic acid
    # Find longest chain starting from acid's alpha carbon
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    chain_lengths = []
    for match in acid_matches:
        alpha_carbon = None
        for neighbor in mol.GetAtomWithIdx(match[0]).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != match[1]:
                alpha_carbon = neighbor.GetIdx()
                break
        if alpha_carbon is None:
            continue
        
        visited = {match[0], alpha_carbon}
        current = alpha_carbon
        length = 1
        
        while True:
            next_carbons = []
            for neighbor in mol.GetAtomWithIdx(current).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    next_carbons.append(neighbor.GetIdx())
            if len(next_carbons) != 1:
                break
            current = next_carbons[0]
            visited.add(current)
            length += 1
        
        chain_lengths.append(length)
    
    if chain_lengths and max(chain_lengths) >= 18:  # Account for possible oxidation shortening
        return True, "C20-derived structure with oxygenation, not prostanoid"
    
    return False, "Doesn't meet C20 backbone requirements"