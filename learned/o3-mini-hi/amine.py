"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms 
by hydrocarbyl groups. Free amines (primary, secondary or tertiary) are those not directly bonded 
to a carbonyl (as in an amide) nor part of a heteroaromatic ring.
 
We improve upon the earlier approach by:
 • More carefully checking if a nitrogen is directly attached to a carbonyl carbon.
 • Skipping only nitrogens that are truly part of a heteroaromatic system (e.g. pyridine‐like, 
   where the nitrogen is sp2 and aromatic) while still allowing, for example, aniline.
 • Using the total hydrogen count from RDKit (GetTotalNumHs) to distinguish primary (–NH2),
   secondary (–NH–) and tertiary amines (–N–).
 
Note: This heuristic does not guarantee perfect classification over all edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_amine(smiles: str):
    """
    Determines if a molecule contains at least one qualifying “free” amine group based on its SMILES.
    A free amine is defined as a nitrogen (atomic number 7) that:
      - Has no positive formal charge.
      - Is not directly bound to a carbonyl carbon (i.e. part of an amide).
      - Is not part of a heteroaromatic ring (e.g. a pyridine‐type N); note that amino groups attached
        to an aromatic ring (e.g. aniline) are allowed if the lone pair remains localized.
      - Has a hydrogen count (from GetTotalNumHs) of 2 (primary), 1 (secondary) or 0 (tertiary).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if at least one qualifying free-amine is found, False otherwise.
        str: Explanation of what was found or why not.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all nitrogen atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # not a nitrogen

        # Reject if positively charged (e.g. quaternary ammonium).
        if atom.GetFormalCharge() > 0:
            continue

        # If the nitrogen is sp2 hybridized and aromatic, it is likely part of a heteroaromatic ring 
        # (e.g. pyridine). Note: aniline-type amines are typically sp3 so they are not excluded.
        if atom.GetHybridization() == Chem.HybridizationType.SP2 and atom.GetIsAromatic():
            continue
        
        # Check if the N is directly bonded to a carbonyl group.
        # That is, for every bond between the N and a carbon, see if that carbon has a double bond to oxygen.
        is_bound_to_carbonyl = False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if other.GetAtomicNum() == 6:  # a carbon neighbor
                # Look over bonds from the carbon to check for a C=O double bond.
                for bond2 in other.GetBonds():
                    # Only consider bonds that are double bonds.
                    if bond2.GetBondType() == rdchem.BondType.DOUBLE:
                        neighbor2 = bond2.GetOtherAtom(other)
                        if neighbor2.GetAtomicNum() == 8:  # oxygen found in a double bond
                            is_bound_to_carbonyl = True
                            break
                if is_bound_to_carbonyl:
                    break
        if is_bound_to_carbonyl:
            # Reject this nitrogen as it is part of an amide (or similar carbonyl linkage).
            continue

        # Determine the number of (implicit+explicit) hydrogen atoms attached.
        nH = atom.GetTotalNumHs()
        if nH == 2:
            return True, "Found a primary amine (-NH2) group not bound to a carbonyl."
        elif nH == 1:
            return True, "Found a secondary amine (-NH-) group not bound to a carbonyl."
        elif nH == 0:
            return True, "Found a tertiary amine (-N-) group not bound to a carbonyl."
    
    return False, "No qualifying free amine functional group found."

# Example usage:
# result, reason = is_amine("CNc1ccccc1")
# print(result, reason)