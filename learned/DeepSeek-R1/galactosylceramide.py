"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:28700 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide (cerebroside with galactose).
    Must have a galactose moiety linked via glycosidic bond to a ceramide structure.
    Ceramide consists of sphingosine backbone with amide-linked fatty acid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Match galactose (beta or alpha, allowing substitutions like sulfates)
    # Pattern captures pyranose ring with at least 4 hydroxyl groups (positions 2,3,4,6)
    galactose = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O")
    if not mol.HasSubstructMatch(galactose):
        return False, "No galactose moiety detected"

    # Find oxygen atom connecting galactose to ceramide (glycosidic bond)
    # Iterate through galactose atoms to find connecting O
    gal_matches = mol.GetSubstructMatches(galactose)
    for match in gal_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0:  # Anomeric carbon
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 2:
                        # Follow the glycosidic bond to ceramide part
                        next_atom = [n for n in neighbor.GetNeighbors() if n.GetIdx() not in match]
                        if not next_atom:
                            continue
                        next_atom = next_atom[0]
                        
                        # Check if connected to sphingosine backbone
                        # Sphingosine pattern: adjacent amine and hydroxyl groups
                        sphingo_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O)-[CH2]-[CH](N)")
                        sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
                        if not sphingo_matches:
                            continue
                        
                        # Verify amide-linked fatty acid
                        # Amide group connected to long aliphatic chain
                        amide_pattern = Chem.MolFromSmarts("[CX4][CX4][NH1][C](=O)[CX4][CX4,H2]")
                        amide_matches = mol.GetSubstructMatches(amide_pattern)
                        if not amide_matches:
                            continue
                        
                        # Check chain lengths (at least 10 carbons in amide chain)
                        for amide in amide_matches:
                            fatty_acid = mol.GetAtomWithIdx(amide[3])  # Carbon in amide
                            chain = AllChem.FindAtomEnvironmentOfRadiusN(mol, 12, fatty_acid.GetIdx())
                            if len(chain) >= 10:  # Approximate chain length check
                                return True, "Galactose linked to ceramide via glycosidic bond"
    
    return False, "No galactose-ceramide structure found"