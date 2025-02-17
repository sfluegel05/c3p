"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol 
Definition: A chromanol with a chroman-6-ol skeleton that is substituted at position 2 
            by a saturated or triply-unsaturated hydrocarbon chain consisting of 
            three isoprenoid units.
            
This program uses rdkit to:
  1. Check that the molecule has a benzopyran (chromanol) core with an aromatic –OH.
  2. Locate the substituent at the “2‐position” (using an atom-mapped SMARTS).
  3. Extract that side chain (ensuring we do not go back into the core) and check that:
       • It is an aliphatic chain composed of carbon only.
       • It is long enough (heuristically, at least 10 carbons).
       • It has either no double bonds (saturated) or exactly three double bonds (triply unsaturated).
       
If any of these checks fail, the function returns False with an explanation.
Note: Because “tocol” chemistry is complex, this is one possible heuristic implementation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_tocol(smiles: str):
    """
    Determines if a molecule (as SMILES string) is a tocol (a vitamin E-like compound)
    based on its chromanol core and a side chain at position 2.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule meets the tocol definition.
       str: A reason explaining the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Step 1. Look for the chromanol (benzopyran) core.
    # Here we use a SMARTS with atom mapping to indicate the substitution at position 2.
    # The pattern "O1CC([*:1])c2ccccc2C1" tries to capture a 2H-1-benzopyran ring system;
    # the atom with map number 1 ([*:1]) is the site for the side chain expected at the 2-position.
    core_smarts = "O1CC([*:1])c2ccccc2C1"
    core_query = Chem.MolFromSmarts(core_smarts)
    core_match = mol.GetSubstructureMatch(core_query)
    
    if not core_match:
        return False, "Chromanol core (benzopyran) not found"

    # Identify the set of atom indices in the core match.
    core_atom_indices = set(core_match)
    
    # Step 2. Check that the benzene (aromatic) part of the core bears an –OH.
    # We use a simple SMARTS to detect an aromatic carbon bearing an OH (phenol).
    phenol_smarts = "c[OH]"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(phenol_smarts)):
        return False, "Aromatic hydroxyl group (phenolic OH) not found in the chromanol core"

    # Step 3. From the mapped core, locate the attachment (side chain) at the 2‑position.
    # We look for the atom in the molecule that was matched to the [*:1] in the SMARTS.
    substitution_root = None
    for atom in mol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            try:
                if int(atom.GetProp("molAtomMapNumber")) == 1:
                    substitution_root = atom
                    break
            except ValueError:
                continue
    if substitution_root is None:
        return False, "Could not locate the 2-position substituent from the core SMARTS match"
        
    # However the substituent root may be “inside” the core query.
    # We want the neighbor of the core that is not part of the core.
    chain_root = None
    for neighbor in substitution_root.GetNeighbors():
        if neighbor.GetIdx() not in core_atom_indices:
            chain_root = neighbor
            break
    if chain_root is None:
        return False, "No side chain found attached at the 2-position"

    # Step 4. Traverse outwards from chain_root (avoiding core atoms) to extract the entire side chain.
    side_chain_idxs = set()
    to_visit = [chain_root.GetIdx()]
    while to_visit:
        current = to_visit.pop()
        if current in side_chain_idxs:
            continue
        side_chain_idxs.add(current)
        current_atom = mol.GetAtomWithIdx(current)
        for nbr in current_atom.GetNeighbors():
            # Skip atoms that are part of the core scaffold.
            if nbr.GetIdx() in core_atom_indices:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())
                
    # Step 5. Check that the side chain is aliphatic; here we insist that all atoms are carbons.
    for idx in side_chain_idxs:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            return False, f"Side chain contains non-carbon atom: {atom.GetSymbol()}"
    
    # Count the number of carbon atoms in the side chain.
    n_chain_carbons = len(side_chain_idxs)
    if n_chain_carbons < 10:
        return False, f"Side chain is too short (found {n_chain_carbons} carbons; expected at least 10 for three isoprenoid units)"
        
    # Step 6. Count the number of double bonds in the side chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        # Consider only bonds entirely within the side chain.
        if a1 in side_chain_idxs and a2 in side_chain_idxs:
            if bond.GetBondType() == BondType.DOUBLE:
                double_bond_count += 1

    # The tocol definition allows either a saturated chain (0 double bonds) or a chain with three C=C.
    if double_bond_count not in (0, 3):
        return False, f"Side chain has {double_bond_count} double bond(s); expected 0 (saturated) or 3 (triply unsaturated)"

    return True, "Molecule contains chromanol core with appropriate aromatic –OH and a 2‑substituent chain matching tocol definition"

# For testing (you can remove these lines when integrating the code):
if __name__ == "__main__":
    test_smiles = "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1"  # (S,R,S)-alpha-tocopherol
    flag, reason = is_tocol(test_smiles)
    print(flag, reason)