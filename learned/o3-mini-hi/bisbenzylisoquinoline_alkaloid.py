"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
#!/usr/bin/env python
"""
Classifies: bisbenzylisoquinoline alkaloid

A bisbenzylisoquinoline alkaloid is defined as a benzylisoquinoline alkaloid whose structures are 
built up of two benzylisoquinoline units linked by ether bridges. Various additional bridging patterns 
(direct carbon–carbon bonds or methylenedioxy groups) can be observed.

Our revised criteria are:
  - The molecule must be a valid SMILES and have molecular weight >= 500 Da.
  - It must contain at least two distinct isoquinoline-like substructures. To capture these,
    we use two SMARTS patterns (ignoring stereochemistry):
      (a) An aromatic isoquinoline: "c1ccc2c(c1)cccn2"
      (b) A tetrahydroisoquinoline-like unit: "c1ccc2CCNC2c1"
    Unique matches from either pattern are collected.
  - There must be a “bridging” motif connecting aromatic portions of two different isoquinoline units.
    We check two possibilities:
      • An ether bridge: an oxygen atom with at least two aromatic neighbors that belong to 
        different isoquinoline units.
      • A methylenedioxy bridge: a CH2 group (sp3 carbon with two neighbors that are oxygen) 
        where each oxygen is connected to an aromatic atom from different isoquinoline units.
  - Finally, the molecule should contain at least 2 nitrogen atoms.
  
If all criteria are met, the function returns (True, <reason>), otherwise (False, <reason>).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    The function checks:
      - Valid SMILES and molecular weight >= 500 Da.
      - Presence of at least two unique isoquinoline-like substructures using two SMARTS patterns.
      - Presence of a bridging motif connecting aromatic substructures belonging to two distinct
        isoquinoline units. Two bridging modes are examined: an ether (oxygen) bridge and a 
        methylenedioxy (-OCH2O-) bridge.
      - At least 2 nitrogen atoms are present.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple with a Boolean decision and an explanation.
    """
    # Convert SMILES to molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to help substructure matching.
    Chem.RemoveStereochemistry(mol)
    
    # Check molecular weight (must be >= 500 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a bisbenzylisoquinoline alkaloid"
    
    # Ensure there are at least 2 nitrogen atoms.
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Not enough nitrogen atoms ({n_count}); at least 2 are needed"
    
    # Define SMARTS patterns (ignoring chirality) for isoquinoline-like substructures.
    arom_iso_smarts = "c1ccc2c(c1)cccn2"  # aromatic isoquinoline
    tetra_iso_smarts = "c1ccc2CCNC2c1"     # tetrahydroisoquinoline-like unit
    arom_iso = Chem.MolFromSmarts(arom_iso_smarts)
    tetra_iso = Chem.MolFromSmarts(tetra_iso_smarts)
    if arom_iso is None or tetra_iso is None:
        return False, "Error creating isoquinoline SMARTS patterns"
    
    # Get substructure matches ignoring chirality.
    arom_matches = mol.GetSubstructMatches(arom_iso, useChirality=False)
    tetra_matches = mol.GetSubstructMatches(tetra_iso, useChirality=False)
    
    # Combine unique substructure matches; each match is a tuple of atom indices.
    unique_iso_units = []
    for match in arom_matches:
        unit = set(match)
        # Avoid duplicates (if the set is already contained in one of the found units)
        if not any(unit <= other for other in unique_iso_units):
            unique_iso_units.append(unit)
    for match in tetra_matches:
        unit = set(match)
        if not any(unit <= other for other in unique_iso_units):
            unique_iso_units.append(unit)
    
    if len(unique_iso_units) < 2:
        return False, f"Found only {len(unique_iso_units)} isoquinoline-like substructure(s); at least 2 are required"
    
    # Define a helper function that, given an atom index, returns the indices of isoquinoline units (by index in unique_iso_units)
    # that contain that atom.
    def iso_units_containing(atom_idx):
        found = set()
        for i, unit in enumerate(unique_iso_units):
            if atom_idx in unit:
                found.add(i)
        return found
    
    # Now search for a valid bridging motif.
    bridge_found = False
    
    # 1. Look for an ether bridge:
    #    Look at every oxygen atom; if it has at least two aromatic neighbors belonging to two different isoquinoline units.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            nei_ids = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIsAromatic()]
            if len(nei_ids) >= 2:
                unit_ids = set()
                for nbr in nei_ids:
                    unit_ids.update(iso_units_containing(nbr))
                if len(unit_ids) >= 2:
                    bridge_found = True
                    break  # one valid bridging motif is enough
    # 2. If not found, look for a methylenedioxy (-OCH2O-) bridge:
    if not bridge_found:
        for atom in mol.GetAtoms():
            # candidate carbon for CH2 bridging: sp3 carbon with atomic number 6,
            # and ideally with 2 neighbors that are oxygen.
            if atom.GetAtomicNum() == 6 and atom.GetHybridization().name == "SP3":
                o_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
                if len(o_neighbors) == 2:
                    # For each oxygen, find aromatic neighbors (other than the bridging carbon)
                    unit_ids = set()
                    for o in o_neighbors:
                        for n in o.GetNeighbors():
                            if n.GetIdx() == atom.GetIdx():
                                continue
                            if n.GetIsAromatic():
                                unit_ids.update(iso_units_containing(n.GetIdx()))
                    if len(unit_ids) >= 2:
                        bridge_found = True
                        break

    if not bridge_found:
        return False, "No valid bridging motif found connecting two distinct isoquinoline units"
    
    # If all criteria pass, return True with a comprehensive explanation.
    reason = ("Molecule has molecular weight {:.1f} Da, contains at least 2 benzylisoquinoline-like substructures "
              "with indices {}, has at least 2 nitrogen atoms, and a valid bridging motif was detected"
              .format(mol_wt, [list(u) for u in unique_iso_units]))
    return True, reason

# For example testing, one can call:
# test_smiles = "[H][C@@]12Cc3ccc(O)c(c3)-c3cc(C[C@]4([H])N(C)CCc5cc(OC)c(OC)c(Oc6cc1c(CCN2C)cc6OC)c45)ccc3OC"  # rodiasine
# print(is_bisbenzylisoquinoline_alkaloid(test_smiles))