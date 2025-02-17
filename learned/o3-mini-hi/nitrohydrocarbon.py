"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C‐nitro compound that is a hydrocarbon in which one or more of the hydrogens 
has been replaced by nitro groups.
In our approach the molecule is classified as a nitrohydrocarbon if:
  1. It contains one or more nitro groups (SMARTS "[N+](=O)[O-]"),
  2. Each nitro group is attached to a carbon atom (i.e. the nitro’s N has at least one C neighbor),
  3. When the nitro groups are ignored, every heavy atom in the scaffold is carbon,
  4. At least one carbon in the scaffold would originally have had a hydrogen attached.
This should rule out species where extra heteroatoms (or excessive substitution) lead to non‐hydrocarbon scaffolds.
"""

from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is defined as a molecule whose core (scaffold), if one were to remove all
    nitro functionalities (i.e. [N+](=O)[O-]), is a pure hydrocarbon (only C and H) and at least one 
    hydrogen remains on the scaffold. In addition, the nitro groups must each be attached to at least
    one carbon atom (i.e. derived from a substitution reaction on a hydrocarbon).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nitrohydrocarbon, False otherwise.
        str: Explanation of the classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find nitro groups using a SMARTS pattern.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Mark all atoms that are part of a nitro group.
    nitro_atom_indices = set()
    valid_nitro_found = False
    for match in nitro_matches:
        # (The match gives indices for atoms in the SMARTS pattern.)
        # Identify the nitrogen atom in the match (atomic number 7).
        nitro_n = None
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                nitro_n = atom
                break
        if nitro_n is None:
            continue
        # For a valid nitro group, the nitrogen should be bonded to at least two oxygens and one carbon.
        neighbors = list(nitro_n.GetNeighbors())
        symbols = [nbr.GetSymbol() for nbr in neighbors]
        if symbols.count("O") >= 2 and symbols.count("C") >= 1:
            valid_nitro_found = True
            nitro_atom_indices.update(match)
    
    if not valid_nitro_found:
        return False, "No nitro group found attached to a carbon atom"
    
    # Check the hydrocarbon scaffold.
    # The scaffold consists of atoms that are not part of any nitro group.
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        # The scaffold must consist only of carbon atoms.
        if atom.GetAtomicNum() != 6:
            return False, f"Found non-hydrocarbon atom ({atom.GetSymbol()}) in the scaffold"
    
    # Ensure that at least one scaffold carbon retains a hydrogen.
    # (i.e. at least one carbon atom not in a nitro group must have at least one attached hydrogen).
    scaffold_has_H = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atom_indices:
            continue
        # GetTotalNumHs counts implicit+explicit hydrogens.
        if atom.GetTotalNumHs() > 0:
            scaffold_has_H = True
            break
    if not scaffold_has_H:
        return False, "No hydrogen remains on the hydrocarbon scaffold"
    
    return True, "Molecule is a nitrohydrocarbon: contains nitro groups attached to a pure hydrocarbon scaffold with remaining hydrogen(s)"

# Example usage:
if __name__ == "__main__":
    examples = {
        "3,7-Dinitrofluoranthene": "[O-][N+](=O)c1cccc-2c1-c1cccc3c(ccc-2c13)[N+]([O-])=O",
        "1-nitroheptane": "[O-][N+](=O)CCCCCCC",
        "Nitrobenzene": "[O-][N+](=O)c1ccccc1",
        "4-nitrotoluene (expected False)": "Cc1ccc(cc1)[N+]([O-])=O",  # per outcome, this example should be rejected
        "Tetranitromethane (expected False)": "C([N+]([O-])=O)([N+]([O-])=O)[N+]([O-])=O"
    }
    for name, smi in examples.items():
        result, reason = is_nitrohydrocarbon(smi)
        print(f"{name}: {result}. Reason: {reason}")