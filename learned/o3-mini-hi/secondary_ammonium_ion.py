"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: Secondary ammonium ion
Definition:
    "An organic cation obtained by protonation of any secondary amino compound; 
     major species at pH 7.3."
A secondary ammonium ion arises from protonation of a secondary amine (R₂NH) to give R₂NH₂⁺.
Thus the key chemical features we use are:
  - A nitrogen (atomic number 7) carrying a formal charge of +1.
  - The nitrogen is sp3-hybridized.
  - It is bonded to exactly 2 heavy (non-hydrogen) atoms.
  - Its total hydrogen count (explicit + implicit) is exactly 2 (so that the sum of bonds equals 4).
  - The overall molecule must be a cation (net charge > 0) so as to not mis‐classify zwitterions.
  - We also screen out cases where a heavy neighbor is “amidic” (i.e. bonded to a double-bonded oxygen) 
    because protonation of an amide is extremely uncommon.
"""

from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if the molecule (given by a SMILES string) contains at least one
    secondary ammonium ion center. Specifically, it looks for a positively charged nitrogen
    (formal charge +1) that is sp3 hybridized and bonded to exactly two heavy (non-hydrogen)
    atoms with a total of 2 hydrogen atoms attached (i.e. mimicking the structure R2NH2+ from a
    protonated secondary amine).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a secondary ammonium ion center is found, False otherwise.
        str: Reason message detailing the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)
    
    # Calculate overall molecule charge by summing formal charge on each atom.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge <= 0:
        return False, "Molecule overall is not a cation (net charge <= 0)."
    
    # Iterate over atoms to find candidate nitrogen centers.
    for atom in mol.GetAtoms():
        # Look for nitrogen with a positive formal charge.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != +1:
            continue
        # Require sp3 hybridization.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue

        # Count total hydrogens attached (explicit+implicit).
        total_H = atom.GetTotalNumHs()
        if total_H != 2:
            continue

        # Now check the environment of each heavy neighbor.
        # We want to avoid cases where a heavy neighbor is part of an amide (i.e. a carbonyl group)
        is_amidic = False
        for nbr in heavy_neighbors:
            # if neighbor is carbon, examine bonds for a double bond to oxygen
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    # Check if the bond is a double bond and the other atom is oxygen.
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_amidic = True
                        break
            if is_amidic:
                break
        if is_amidic:
            continue  # Skip candidate if one substituent appears carbonyl-like.
        
        # If we reach here, we have a candidate center meeting our criteria.
        return True, ("Found a sp3 nitrogen with formal charge +1, bonded to exactly 2 heavy atoms "
                       "and 2 hydrogens. This is consistent with a protonated secondary amine (secondary ammonium ion).")
    
    return False, "No secondary ammonium ion center found based on the defined criteria."

# If this script is run as main, include a few tests.
if __name__ == "__main__":
    test_smiles = [
        # True positives (some examples from the provided list)
        "[H]\\C(C)=C1/C[NH2+][C@@]2([H])CC3=C(NC4=CC=CC=C34)C(=O)C[C@]1([H])[C@]2([H])C(=O)OC",  # perivine(1+)
        "C1COCC[NH2+]1",  # morpholinium
        "C[NH2+]C",  # dimethylaminium
        "C=1C=CC=C2C(=CNC12)CC[NH2+]C",  # N-methyltryptaminium
        # A false positive example: a zwitterion should be rejected.
        "O[C@@H]1C[NH2+][C@@H](C1)C([O-])=O",  # cis-4-hydroxy-L-proline zwitterion
    ]
    for s in test_smiles:
        result, reason = is_secondary_ammonium_ion(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")