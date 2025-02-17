"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: Secondary ammonium ion
Definition:
    "An organic cation obtained by protonation of any secondary amino compound; 
     major species at pH 7.3."

A secondary ammonium ion arises from protonation of a secondary amine (R2NH) to give R2NH2+.
Key chemical features used here are:
  - A nitrogen (atomic number 7) carrying a formal charge of +1.
  - The nitrogen (after adding hydrogens) must have a total of 4 bonds.
  - Exactly 2 of these bonds are made to heavy (non-hydrogen) atoms.
  - Total number of attached hydrogens (explicit + implicit) is 2.
  - The overall molecule must be a cation (net charge > 0).
  - We screen out cases where a heavy neighbor is “amidic” (i.e. has a double bond to oxygen)
    on a bond other than the one to the candidate nitrogen.
"""

from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if the molecule (given by a SMILES string) contains at least one
    secondary ammonium ion center. Specifically, it looks for a nitrogen atom that is:
      - Positively charged (+1 formal charge)
      - sp3 hybridized
      - Has been protonated so that its total connectivity (after adding hydrogens) is 4,
        with exactly 2 of these bonds going to non-hydrogen (heavy) atoms and 2 hydrogens.
      - Its heavy neighbors (non-hydrogen atoms) do not have a double bond to oxygen 
        (aside from the bond they share with the candidate nitrogen).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a secondary ammonium ion center is found, False otherwise.
        str: Reason message detailing the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to accurately count attached H's.
    mol = Chem.AddHs(mol)
    
    # Calculate overall molecular charge: must be > 0.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge <= 0:
        return False, "Molecule overall is not a cation (net charge <= 0)."
    
    # Loop over atoms to find candidate ammonium centers.
    for atom in mol.GetAtoms():
        # Check if atom is nitrogen with positive charge.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != +1:
            continue
        
        # Ensure nitrogen is sp3 hybridized.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # After adding H's, the candidate ammonium ion must be tetra-coordinated.
        if atom.GetDegree() != 4:
            continue
        
        # Count heavy (non-H) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue

        # Count total hydrogens attached (explicit + implicit)
        total_h = atom.GetTotalNumHs()
        if total_h != 2:
            continue

        # Check each heavy neighbor's environment to screen out "amidic" groups.
        # We want to ignore bonds that are to the candidate nitrogen.
        amidic_flag = False
        for nbr in heavy_neighbors:
            for bond in nbr.GetBonds():
                # Skip the bond that connects this heavy neighbor to our candidate N.
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # If neighbor is carbon, and has a double bond (C=O) with oxygen, mark it as amidic.
                if nbr.GetAtomicNum() == 6:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        amidic_flag = True
                        break
            if amidic_flag:
                break
        if amidic_flag:
            continue

        # If a candidate ammonium ion center is found, return success.
        return True, ("Found a sp3 nitrogen with formal charge +1 that is tetra-coordinated "
                       "with exactly 2 heavy atoms and 2 hydrogens. This is consistent with a "
                       "protonated secondary amine (secondary ammonium ion).")
        
    return False, "No secondary ammonium ion center found based on the defined criteria."

# If this script is run as the main program, run a few tests.
if __name__ == "__main__":
    test_smiles = [
        # Example structures belonging to the secondary ammonium ion class:
        "[H]\\C(C)=C1/C[NH2+][C@@]2([H])CC3=C(NC4=CC=CC=C34)C(=O)C[C@]1([H])[C@]2([H])C(=O)OC",  # perivine(1+)
        "C1COCC[NH2+]1",  # morpholinium
        "CC(C[NH2+]C)C",  # N-methyl-2-methylpropan-1-aminium
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine(1+)
        "[H][C@]12CN(C)[C@]([H])(C[NH2+]1)CC1=C[C@]([H])(C(=O)CC1)[C@]1([H])C=C(CCC1=O)C2",  # herquline B(1+)
        "C1=C(C(=CC2=C1[C@H]([NH2+]CC2)CC3=CC=C(C(=C3)O)OC)OC)O",  # (R)-norreticuline(1+)
        "C=1C=CC=C2C(=CNC12)CC[NH2+]C",  # N-methyltryptaminium
        "C[NH2+]C",  # dimethylaminium
        "[H][C@]12C[C@@]3([H])[C@]4([H])CC=C5C[C@H](CC[C@]5(C)[C@@]4([H])CC[C@]3(C)[C@@]1([H])[C@H](C)[C@@]1(CC[C@@H](C)C[NH2+]1)O2)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",  # solasodine 3-beta-D-glucoside(1+)
        "[C@]12(C(CCC(=C1)C[C@H]3CN[C@H](C[NH2+]3)CC4=C[C@]2(C(CC4)=O)[H])=O)[H]",  # (1S,7R,8R,14S)-6,9-dioxo-15,17-diazatetracyclo[...] (truncated)
        "[H][C@]12CN(C)[C@]([H])(C[NH2+]1)CC1=C[C@@]([H])(C(=O)CC1)[C@@]1([H])C=C(CCC1=O)C2",  # herquline C(1+)
        "C=1C=C(OCC[NH2+][C@H](CC2=CC=C(C(=C2)S(N)(=O)=O)OC)C)C(=CC1)OCC",  # ent-tamsulosin(1+)
        "C1C[C@]2([C@@](CC[C@@]3([C@@]2(CC[C@@H](C3)O)C)[H])([C@]4([C@]1([C@@]5([C@](C4)(O[C@@]6([C@H]5C)[NH2+]C[C@@H](C)CC6)[H])[H])C)[H])[H])[H]",  # tomatidine(1+)
        "OC1=CC2=C(C[NH2+]CC2)C=C1O",  # norsalsolinol(1+)
        "C[NH2+]CC1=CC=C(C=C1)C1=CC2=CC(F)=CC(C(N)=O)=C2O1",  # mefuparib(1+)
        "C[NH2+]CCC1=CC(OC)=C(OC)C(OC)=C1",  # N-methylmescalinium
        "C[C@H](CCC1=CC=CC=C1)[NH2+]C[C@H](O)C1=CC(C(N)=O)=C(O)C=C1",  # dilevalol(1+)
        "C(C[NH3+])CCC[NH2+]CC([C@@H](CO)O)=O",  # N-(D-erythrulosyl)-cadaverine(2+)
        "[C@@H]12[C@H](O[C@@]3(O[C@H](C)CC([C@@]3(O1)O)=O)[H])[C@H]([C@H]([NH2+]C)[C@@H]([C@@H]2NC)O)O",  # spectinomycin(1+)
        "CCC[C@H]([NH2+][C@H]1CCC2=CC(F)=CC(F)=C2C1)C(=O)NC1=CN(C=N1)C(C)(C)C[NH2+]CC(C)(C)C",  # nirogacestat(2+)
        "CC(C)[NH2+]CC(Cl)C1=CC2=CC=CC=C2C=C1",  # N-[2-chloro-2-(naphthalen-2-yl)ethyl]propan-2-aminium
        "C(CC[NH3+])CC[NH+](C)C",  # N,N-dimethyl-cadaverine(+2)
        "C1=C(C(=CC2=C1C([NH2+]CC2)CC3=CC=C(C(=C3)O)OC)OC)O",  # norreticuline(1+)
        "O=S(=O)(N1CCC[NH2+]CC1)C1=CC=CC2=C1C=CN=C2",  # fasudil(1+)
        "CC[NH2+]CC",  # diethylammonium
        "C1=C(C(=CC2=C1[C@@H]([NH2+]CC2)CC3=CC=C(C(=C3)O)OC)OC)O",  # (S)-norreticuline(1+)
        "C1CCCC[NH2+]1",  # piperidinium
    ]
    
    for smiles in test_smiles:
        result, reason = is_secondary_ammonium_ion(smiles)
        print("SMILES: {}\nResult: {}\nReason: {}\n".format(smiles, result, reason))