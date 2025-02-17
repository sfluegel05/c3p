"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: Mineral (inorganic, often ionic/geologically formed species)
Based on heuristics:
 - Most minerals occur as disconnected fragments (ionic salts)
 - Many carry explicit charges.
 - They often include elements that are not typical for organic substances.
Examples include calcium difluoride, disodium tetraborate, pyrite, etc.
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is considered a mineral based on its SMILES string.
    Heuristic criteria:
      1. The substance must be parseable.
      2. Many minerals are salts or ionic: if the molecule consists of multiple disconnected fragments, this suggests ionic mineral formation.
      3. Many minerals contain ions, i.e. one or more atoms carry non-zero formal charges.
      4. Inorganic substances typically include elements beyond the limited set common to organic molecules.
         We define the common organic elements as: H, C, N, O, F, Cl, Br, I, P, and S.
         Any atom outside this set (e.g. metals or metalloids like Ca, Mg, Fe, Si, Al, Sb, etc.) increases the odds that the substance is mineral.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fulfills our criteria for a mineral, False otherwise.
        str: A reason describing the classification.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for multiple fragments (common in salts/ionic minerals)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return True, f"Multiple fragments detected (count: {len(frags)}), typical for ionic mineral compounds."
    
    # Check if any atom has a non-zero formal charge (suggesting an ionic component)
    if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return True, "Presence of ions with non-zero formal charges suggests a mineral (ionic compound)."
    
    # Define a set of elements common in organic compounds.
    # (Atomic numbers: H=1, C=6, N=7, O=8, F=9, P=15, S=16, Cl=17, Br=35, I=53)
    organic_set = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    
    # Examine each atom: if any atom is not in the organic set, it likely indicates an inorganic (mineral) nature.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in organic_set:
            return True, f"Contains element {atom.GetSymbol()} (atomic number {atom.GetAtomicNum()}) not typical in organic compounds."
    
    # If none of the above conditions are met, then the molecule is likely a simple covalent (organic-like) species.
    return False, "Does not display inorganic or ionic features typical of minerals."

# Example calls (you may uncomment and test with provided SMILES examples):
# print(is_mineral("[S--].[S--].[S--].[S--].[Fe++].[Fe+3].[Fe+3]"))  # greigite
# print(is_mineral("[F-].[F-].[Ca++]"))  # calcium difluoride
# print(is_mineral("Cl[Cu]Cl.O.O.O.O.O"))  # copper(II) chloride pentahydrate
# print(is_mineral("[Sb](F)(F)(F)(F)F"))  # antimony pentafluoride