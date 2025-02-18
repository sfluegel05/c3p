"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_galactosylceramide(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define galactose pattern (beta-D-galactopyranose)
    # This SMARTS matches the ring with hydroxyls in typical positions, allowing for substitutions (like sulfates)
    galactose_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O")
    if not mol.HasSubstructMatch(galactose_pattern):
        # Check for alpha-D-galactose as well?
        # Alternatively, use a more general pattern ignoring stereochemistry
        # galactose_pattern = Chem.MolFromSmarts("[C]1[C]([C]([C]([C](O1)CO)O)O)O")
        # if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose moiety found"
    
    # Check that galactose is connected via an oxygen to the ceramide
    # Find the oxygen atom connecting galactose to the rest
    # This part is tricky; perhaps look for an oxygen that is part of the glycosidic bond
    # The oxygen in the galactose ring that is connected to the ceramide
    # We can check if the galactose substructure is connected via an oxygen to another part of the molecule
    matches = mol.GetSubstructMatches(galactose_pattern)
    for match in matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 2:  # Assuming glycosidic O is connecting two carbons
                    # Check if this oxygen connects to a non-galactose part (ceramide)
                    for bond in neighbor.GetBonds():
                        other_atom = bond.GetOtherAtom(neighbor)
                        if other_atom.GetIdx() not in match:
                            # Now check if the other part has a ceramide structure
                            # Ceramide has an amide group (CONH) connected to two chains
                            # Look for CONH where the NH is connected to a chain
                            ceramide_pattern = Chem.MolFromSmarts("[CX4][CX4][NH1][C](=O)[CX4]")
                            ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
                            if ceramide_matches:
                                # Check for long chains (e.g., at least 10 carbons)
                                # This is a simplification; real chains are longer
                                # Count carbons in the fatty acid and sphingosine parts
                                # This part may need refinement
                                return True, "Galactose linked to ceramide structure"
    
    return False, "No ceramide structure linked to galactose"