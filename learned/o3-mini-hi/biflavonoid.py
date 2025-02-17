"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: Biflavonoid

Definition: A flavonoid oligomer obtained by the oxidative coupling of at least 
two units of aryl‐substituted benzopyran rings (or their substituted derivatives),
resulting in the two ring systems being joined together by a single atom or bond.

This improved heuristic performs the following:
  • Basic sanity checks: the molecule must have a minimum number of rings (≥6) and a 
    molecular weight >400 Da.
  • It then scans for single bonds connecting two aromatic atoms.
  • For each candidate bond, the molecule is “broken” into two fragments.
  • Each fragment is tested for the presence of a flavonoid‐like substructure via 
    three SMARTS patterns (flavone‐like, flavan‐like, or flavanone‐like cores) plus a 
    requirement that the fragment contains at least two rings.
  • If a bond is found for which both fragments qualify, the molecule is classified as a biflavonoid.
Note: This remains a heuristic classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is considered a biflavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a biflavonoid, else False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Global sanity: require a minimum number of rings and sufficient molecular weight.
    num_rings = len(mol.GetRingInfo().AtomRings())
    if num_rings < 6:
        return False, f"Too few rings ({num_rings}) to be a biflavonoid"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a biflavonoid"
    
    # Define SMARTS patterns for flavonoid-like cores.
    # Three patterns are used:
    #  1. Flavone-like core (chromen-4-one)
    #  2. Flavan-like core (saturated at C2 of the pyran ring)
    #  3. Flavanone-like core (phenyl attached to a carbonyl-bearing pyran ring)
    smarts_list = [
        "c1ccc2c(c1)oc(=O)c(c2)",    # flavone-like unit
        "c1ccc2c(c1)OC(c2)",          # flavan-like unit
        "c1ccc(C2=CC(=O)OC2)cc1"       # flavanone-like unit
    ]
    patterns = []
    for sp in smarts_list:
        patt = Chem.MolFromSmarts(sp)
        if patt is None:
            return False, f"Error generating pattern from SMARTS: {sp}"
        patterns.append(patt)
    
    def has_flavonoid_unit(fragment):
        """
        Checks if a given fragment contains at least one flavonoid-like core.
        In addition to matching one of the SMARTS patterns, we require that the fragment has at least 2 rings.
        """
        # Check for minimal ring count in the fragment (flavonoid cores generally have 2–3 rings)
        if rdMolDescriptors.CalcNumRings(fragment) < 2:
            return False
        for patt in patterns:
            if fragment.HasSubstructMatch(patt):
                return True
        return False
    
    # Scan for bonds that may link two flavonoid units.
    # We look for a single bond that connects two aromatic atoms.
    for bond in mol.GetBonds():
        # Only consider single bonds
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        # Only consider bonds between two aromatic atoms
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if not (atom1.GetIsAromatic() and atom2.GetIsAromatic()):
            continue
        
        # Create a copy of the molecule as a mutable RWMol and remove the candidate bond.
        rwmol = Chem.RWMol(mol)
        rwmol.RemoveBond(atom1.GetIdx(), atom2.GetIdx())
        # Get the disconnected fragments as separate molecules
        frags = Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=True)
        if len(frags) < 2:
            continue
        
        # Check if at least two fragments have a flavonoid-like substructure.
        flav_count = 0
        for frag in frags:
            if has_flavonoid_unit(frag):
                flav_count += 1
        if flav_count >= 2:
            return True, f"Found a linking bond connecting flavonoid units; molecular weight = {mol_wt:.1f} Da"
    
    return False, "No linking bond detected that connects two flavonoid-like units"

# Simple testing with an example from the provided list.
if __name__ == "__main__":
    # Example: (+)-morelloflavone (a known biflavonoid)
    test_smiles = "Oc1ccc(cc1)[C@@H]1Oc2cc(O)cc(O)c2C(=O)[C@H]1c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1"
    res, reason = is_biflavonoid(test_smiles)
    print("Biflavonoid:", res)
    print("Reason:", reason)