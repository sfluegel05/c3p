"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
Definition: Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones,
derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups
(with any necessary rearrangement of double bonds). (Polycyclic and heterocyclic analogues are included.)

This improved approach:
  - Uses RDKit’s SSSR (unique rings) to iterate over rings.
  - Only considers rings that are at least 5-membered and are composed entirely of carbon atoms that seem sp2 (or aromatic).
  - For each ring, counts “exocyclic carbonyls”: ring carbon atoms with a C=O bond where the oxygen is not part of the ring.
  - For 6-membered rings (typical quinones), further requires that at least one pair of carbonyls is separated by three bonds along the ring 
    (as in p-benzoquinones). Rings of other sizes only need to have at least two such carbonyl groups.
  - Returns a boolean and an explanatory reason.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.

    Strategy:
      1. Parse the SMILES.
      2. Get the set of SSSR rings (which represent unique rings even in fused systems).
      3. For each ring with at least 5 atoms:
           - Check that every atom in the ring is carbon and has sp2 geometry or is aromatic.
           - Identify exocyclic carbonyl substituents: for a carbon in the ring, check if it is double‐bonded to an oxygen
             that is not also in the ring.
           - For 6-membered rings, further require that at least one pair of carbonyl-containing atoms is separated by three positions
             (i.e. “para” on the ring) to mirror the classical quinone structure.
           - For rings of other sizes, require at least two carbonyls.
      4. Return True, with a formatted reason if any ring meets these criteria; otherwise return False.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a quinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain SSSR rings (each as an ordered tuple of atom indices)
    ssr_rings = list(Chem.GetSymmSSSR(mol))
    if not ssr_rings:
        return False, "No rings found in molecule"

    for ring in ssr_rings:
        # Only consider rings with at least 5 atoms.
        if len(ring) < 5:
            continue

        # Get the ordered list of atom indices for the ring.
        ring_indices = list(ring)
        
        # Check that every atom is carbon and is sp2 or aromatic.
        all_carbon_sp2 = True
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                all_carbon_sp2 = False
                break
            # Check for aromaticity or sp2 hybridization.
            if not (atom.GetIsAromatic() or atom.GetHybridization() == rdchem.HybridizationType.SP2):
                all_carbon_sp2 = False
                break
        if not all_carbon_sp2:
            continue

        # Now look for exocyclic carbonyl groups attached to ring carbons.
        carbonyl_atoms = []  # store indices of ring carbons that bear an external C=O
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            has_carbonyl = False
            for bond in atom.GetBonds():
                # Only consider double bonds.
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                nei = bond.GetOtherAtom(atom)
                # Check if neighbor is oxygen.
                if nei.GetAtomicNum() != 8:
                    continue
                # Exclude case where oxygen is part of the same ring.
                if nei.GetIdx() in ring_indices:
                    continue
                has_carbonyl = True
                break
            if has_carbonyl:
                carbonyl_atoms.append(idx)
                
        if len(carbonyl_atoms) < 2:
            continue
        
        # For 6-membered rings, further check that at least one pair of these carbonyl-bearing atoms
        # is separated by 3 positions in the ring ordering (i.e. para relationship).
        if len(ring_indices) == 6:
            # First, try to get an ordering. (Chem.GetSymmSSSR returns rings in arbitrary order,
            # but for a simple cycle we can assume the indices list is cyclically ordered.)
            found_para = False
            # Create a mapping from ring position to atom index.
            # We assume the list order is consistent; if not, we can try to order by connectivity.
            n = 6
            # Find positions of carbonyl atoms in the ring ordering.
            positions = []
            for pos, idx in enumerate(ring_indices):
                if idx in carbonyl_atoms:
                    positions.append(pos)
            # Check all unique pairs to see if cyclic distance is 3.
            for i in range(len(positions)):
                for j in range(i+1, len(positions)):
                    d = abs(positions[i] - positions[j])
                    cyclic_distance = min(d, n - d)
                    if cyclic_distance == 3:
                        found_para = True
                        break
                if found_para:
                    break
            if not found_para:
                continue  # skip this ring if no proper spacing of carbonyls
        
        # Formulate a reason string.
        reason = (f"Found ring with atoms {tuple(ring_indices)}: {len(carbonyl_atoms)} exocyclic carbonyl groups "
                  f"attached to an all-carbon conjugated ring of size {len(ring_indices)}.")
        return True, reason

    return False, "No fully conjugated all‐carbon ring with at least two exocyclic carbonyl groups (with proper spacing in 6‑membered rings) was found."


# Example calls (uncomment to test)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin (true case)
# print(is_quinone("O=C1C(=C[C@@H]2C=C[C@@H]3CC(N[C@]43[C@@H]2C(=O)C5=C(C(O)=C(C)C(=C5C4=O)O)C([C@@H](C=C[C@H]([C@@H]1C)O)C)=O)=O)CC(C)C"))  # Ansaseomycin A (example false negative in previous attempt)