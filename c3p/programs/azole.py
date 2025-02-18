"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
Defined as any monocyclic heteroarene consisting of a five‐membered aromatic ring 
containing at least one nitrogen. The ring can also contain other non-carbon atoms 
such as O or S. In order to reduce false positives (e.g. peptides or large multi‐ring 
molecules that merely contain an azole fragment), the algorithm also checks that the 
candidate azole ring is not a very minor substructure and that the molecule does 
not contain multiple amide bonds (a hint for a peptide).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    
    An azole is defined here as a molecule that contains at least one five‐membered 
    aromatic heterocycle containing at least one nitrogen, and in which that ring is not merely 
    a small substituent in a larger system. For example, if the candidate azole ring makes up only a very small 
    fraction of the heavy atoms or if the molecule contains multiple amide bonds (as occurs in peptides),
    then the molecule is not classified as an azole.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an azole, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"

    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"

    # Allowed atomic numbers for the azole ring: carbon, nitrogen, oxygen, sulfur.
    allowed_atoms = {6, 7, 8, 16}
    
    candidate_rings = []
    # Loop over each ring found in the molecule.
    for ring in atom_rings:
        # We are interested only in five-membered rings.
        if len(ring) != 5:
            continue
        
        # Get the atoms in the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check that every atom in the ring is aromatic.
        if not all(atom.GetIsAromatic() for atom in ring_atoms):
            continue
        
        # Check that every atom is one of the allowed elements.
        if not all(atom.GetAtomicNum() in allowed_atoms for atom in ring_atoms):
            continue
        
        # Check that there is at least one nitrogen in the ring.
        if not any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            continue
        
        # This ring qualifies as a candidate five-membered heteroaromatic ring.
        candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No five-membered aromatic heterocycle with nitrogen found"

    # Heuristic 1: Check for multiple amide bonds.
    # (Multiple amide bonds often indicate a peptide or polypeptide rather than a single azole moiety.)
    amide_smarts = Chem.MolFromSmarts("[C](=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) > 1:
        return False, f"Molecule contains {len(amide_matches)} amide bonds, likely a peptide"

    # Heuristic 2: Check the relative size of each candidate ring compared
    # to the total number of heavy atoms.
    num_heavy = mol.GetNumHeavyAtoms()
    # We record the maximum fraction among candidate rings.
    max_fraction = 0.0
    for ring in candidate_rings:
        fraction = len(ring) / num_heavy
        if fraction > max_fraction:
            max_fraction = fraction
    
    # If the candidate ring accounts for very little of the molecule, then it may only be a substructure.
    # (Threshold value: 0.15 – this can be tuned. For example, in atorvastatin the five‐membered ring
    # may be a small part of a larger scaffold. In such cases the azole might still be significant. Adjust as needed.)
    if max_fraction < 0.15:
        return False, f"Candidate azole ring accounts for only {max_fraction*100:.1f}% of the molecule's heavy atoms"
    
    # If all checks pass, we classify the molecule as an azole.
    return True, "Found a five‐membered aromatic heterocycle with nitrogen (azole ring)"

# Example usage: (For production you may remove or comment out this section)
if __name__ == "__main__":
    examples = [
        ("CC(C)OCCCNC1=NC=NC2=C1C=NN2C3=CC=CC=C3", "1-phenyl-N-(3-propan-2-yloxypropyl)-4-pyrazolo[3,4-d]pyrimidinamine"),
        ("Cn1cnc(CCN)c1", "N(tele)-methylhistamine"),
        ("O1C(=NC=C1)CCCCC", "2-Pentyloxazole"),
        ("O=C(N[C@@H](CO)C(O)=O)[C@@H](NC(=O)[C@@H](N)C)CC1NC=NC1", "Ala-His-Ser")
    ]
    for smi, name in examples:
        result, reason = is_azole(smi)
        print(f"Test: {name}, Result: {result}, Reason: {reason}")