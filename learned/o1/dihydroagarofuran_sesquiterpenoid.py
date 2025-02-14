"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""

from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid has a fused tricyclic system consisting of two cyclopentane rings
    and a tetrahydrofuran ring (a five-membered ring containing an oxygen atom), forming the dihydroagarofuran skeleton.
    Additionally, it is a sesquiterpenoid (contains 15 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdMolDescriptors
    
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sesquiterpenoid (15-carbon terpenoid)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Number of carbons ({num_carbons}) less than 15"

    # Get ring information
    ssr = Chem.GetSymmSSSR(mol)
    atom_rings = [set(ring) for ring in ssr]

    # Identify five-membered rings containing an oxygen atom (THF rings)
    thf_rings = []
    for ring in atom_rings:
        if len(ring) == 5:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_nums = [atom.GetAtomicNum() for atom in ring_atoms]
            if 8 in atom_nums:  # Oxygen atom present
                thf_rings.append(ring)

    if not thf_rings:
        return False, "No tetrahydrofuran rings found"

    # For each THF ring, check if it's fused to exactly two cyclopentane rings
    for thf_ring in thf_rings:
        fused_rings = []
        for other_ring in atom_rings:
            if other_ring == thf_ring:
                continue
            # Check if rings are fused (share at least two atoms)
            shared_atoms = thf_ring.intersection(other_ring)
            if len(shared_atoms) >= 2:
                fused_rings.append(other_ring)

        # Check if there are exactly two fused rings
        if len(fused_rings) != 2:
            continue

        # Check if both fused rings are five-membered carbocyclic rings (cyclopentane)
        valid_fused = True
        for ring in fused_rings:
            if len(ring) != 5:
                valid_fused = False
                break
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_nums = [atom.GetAtomicNum() for atom in ring_atoms]
            if any(an != 6 for an in atom_nums):  # Not all carbons
                valid_fused = False
                break

        if valid_fused:
            # Found a dihydroagarofuran skeleton
            return True, "Contains dihydroagarofuran skeleton with fused tricyclic system"

    return False, "Dihydroagarofuran skeleton not found"