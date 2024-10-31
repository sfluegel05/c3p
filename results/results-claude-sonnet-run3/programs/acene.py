from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetSymmSSSR

def is_acene(smiles: str):
    """
    Determines if a molecule is an acene (linearly fused benzene rings).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an acene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check molecule is made only of carbons and hydrogens
    atoms = mol.GetAtoms()
    if not all(atom.GetSymbol() in ['C', 'H'] for atom in atoms):
        return False, "Contains non C/H atoms"

    # Get all rings
    rings = list(GetSymmSSSR(mol))
    
    # Check all rings are 6-membered and aromatic
    if not all(len(ring) == 6 for ring in rings):
        return False, "Contains non 6-membered rings"
        
    ring_atoms = set()
    for ring in rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetIsAromatic() for atom in atoms):
            return False, "Contains non-aromatic rings"
        ring_atoms.update(ring)

    # Create graph of ring connectivity
    ring_graph = {i:[] for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if len(set(rings[i]).intersection(set(rings[j]))) == 2:
                ring_graph[i].append(j)
                ring_graph[j].append(i)

    # Check if rings form a linear chain
    # Each ring should connect to max 2 other rings
    if not all(len(connections) <= 2 for connections in ring_graph.values()):
        return False, "Rings not in linear arrangement"

    # Find endpoints (rings connected to only 1 other ring)
    endpoints = [ring for ring in ring_graph if len(ring_graph[ring]) == 1]
    if len(endpoints) != 2:
        return False, "Not a linear chain of rings"

    # Try to traverse from one endpoint to the other
    current = endpoints[0]
    path = [current]
    while ring_graph[current]:
        next_options = [r for r in ring_graph[current] if r not in path]
        if not next_options:
            break
        current = next_options[0]
        path.append(current)

    if len(path) != len(rings):
        return False, "Rings not in linear arrangement"

    # Additional check for non-linear fusion pattern
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            shared_atoms = set(rings[i]).intersection(set(rings[j]))
            if len(shared_atoms) > 0 and j - i > 1:
                return False, "Contains non-linear ring fusion"

    acene_names = {
        1: "benzene",
        2: "naphthalene", 
        3: "anthracene",
        4: "tetracene",
        5: "pentacene",
        6: "hexacene",
        7: "heptacene",
        8: "octacene",
        9: "nonacene",
        10: "decacene"
    }

    name = acene_names.get(len(rings), f"{len(rings)}-acene")
    return True, f"Structure is {name} ({len(rings)} linearly fused benzene rings)"
# Pr=0.6363636363636364
# Recall=1.0