"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms 
by hydrocarbyl groups. Primary (–NH2), secondary (–NH–) or tertiary (–N–) amines are valid.
A nitrogen directly bonded to a carbonyl (e.g. in an amide) or that is a part of a heteroaromatic
ring (e.g. pyridine) is excluded.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_amine(smiles: str):
    """
    Determines if a molecule contains a qualifying amine functional group based on its SMILES string.
    The algorithm iterates over nitrogen atoms and checks:
      - The nitrogen atom is not quaternary.
      - It is not part of a heteroaromatic ring (e.g. pyridine).
      - It is not directly attached to a carbonyl moiety (i.e. part of an amide).
      - It has a hydrogen count consistent with a primary (2 H), secondary (1 H) or tertiary (0 H) 
        derivative of ammonia.
        
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if at least one qualifying free-amine group is found, False otherwise.
        str: Explanation of what was found or why not.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms looking for nitrogen atoms (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # skip non-nitrogen atoms

        # Skip atoms with a positive formal charge (e.g. quaternary ammonium)
        if atom.GetFormalCharge() > 0:
            continue

        # Exclude nitrogen atoms in aromatic rings (e.g. pyridine); they are not considered free amines.
        if atom.GetIsAromatic() and atom.IsInRing():
            continue
        
        # Also reject nitrogens that appear over-coordinated
        if atom.GetDegree() > 3:
            continue
        
        # Check for attachment to a carbonyl group.
        attached_to_carbonyl = False
        for neighbor in atom.GetNeighbors():
            # Only consider carbon neighbors
            if neighbor.GetAtomicNum() == 6:
                # Check each bond of the neighboring carbon.
                for bond in neighbor.GetBonds():
                    # Look for a double bond from the carbon neighbor to an oxygen.
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(neighbor)
                        if other.GetAtomicNum() == 8:
                            attached_to_carbonyl = True
                            break
                if attached_to_carbonyl:
                    break
        if attached_to_carbonyl:
            continue

        # Get the total number of bonded hydrogens on this nitrogen.
        nH = atom.GetTotalNumHs()
        
        # Classify the amine type based on hydrogen count
        if nH == 2:
            return True, "Found a primary amine (-NH2) group not bound to a carbonyl."
        elif nH == 1:
            return True, "Found a secondary amine (-NH-) group not bound to a carbonyl."
        elif nH == 0:
            return True, "Found a tertiary amine (-N-) group not bound to a carbonyl."

    # No qualifying free amine functional group was found.
    return False, "No qualifying free amine functional group found."

# Example usage:
# result, reason = is_amine("CNc1ccccc1")
# print(result, reason)