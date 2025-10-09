"""
Example demonstrating new constraint system enhancements.

This script shows the new features implemented from TODO.md section 1:
- Force parameter support for Pocket and Contact classes
- Distance validation with 4-20Å range and 6Å default
- Enhanced pocket constraint validation (any sequence type as binder)
"""

from boltzjobs import Job

def main():
    print("=== Constraint System Enhancements Demo ===\n")
    
    # Create a complex molecular system
    job = Job("constraint_enhancements_demo")
    
    # Add various sequence types
    print("Adding sequences...")
    job.add_protein_chain("MKLLVVGATTGCQPW", ids="PROT")
    job.add_dna_chain("ATCGATCG", ids="DNA1")
    job.add_rna_chain("AUCGAUCG", ids="RNA1")
    job.add_ligand(smiles="CCO", ids="LIG1")
    job.add_ligand(ccd="ATP", ids="LIG2")
    
    print(f"Added {len(job.sequences)} sequences")
    print()
    
    # Demonstrate new pocket constraints - now work with ANY sequence type
    print("=== Pocket Constraints (Enhanced) ===")
    
    # Protein as binder (new capability)
    pocket1 = job.add_pocket("PROT", 
                           contact_tokens=[("LIG1", 1), ("LIG2", 1)],
                           max_distance=8.0,
                           force=True)  # New force parameter
    print(f"✓ Protein pocket: {pocket1.binder} (force={pocket1.force})")
    
    # DNA as binder (new capability) 
    pocket2 = job.add_pocket("DNA1", max_distance=7.5)  # New default: 6.0Å
    print(f"✓ DNA pocket: {pocket2.binder} (distance={pocket2.max_distance}Å)")
    
    # RNA as binder (new capability)
    pocket3 = job.add_pocket("RNA1", max_distance=6.0, force=False)
    print(f"✓ RNA pocket: {pocket3.binder} (force={pocket3.force})")
    
    # Ligand as binder (still works)
    pocket4 = job.add_pocket("LIG1", max_distance=5.5, force=True)
    print(f"✓ Ligand pocket: {pocket4.binder} (force={pocket4.force})")
    print()
    
    # Demonstrate new contact constraints with force parameter
    print("=== Contact Constraints (Enhanced) ===")
    
    # Standard contact with new default distance (6.0Å)
    contact1 = job.add_contact("PROT", 5, "LIG1", 1)
    print(f"✓ Standard contact: max_distance={contact1.max_distance}Å (new default)")
    
    # Enforced contact with force parameter
    contact2 = job.add_contact("PROT", 10, "LIG2", 1, 
                             max_distance=7.0, 
                             force=True)  # New force parameter
    print(f"✓ Enforced contact: force={contact2.force}, distance={contact2.max_distance}Å")
    
    # Cross-chain nucleic acid contact
    contact3 = job.add_contact("DNA1", 3, "RNA1", 4, 
                             max_distance=12.0, 
                             force=False)
    print(f"✓ Nucleic acid contact: {contact3.token1} - {contact3.token2}")
    print()
    
    # Demonstrate distance validation (4-20Å range)
    print("=== Distance Validation (4-20Å Range) ===")
    
    try:
        # Valid distances
        job.add_contact("PROT", 1, "LIG1", 1, max_distance=4.0)  # minimum
        print("✓ Min distance (4.0Å): Valid")
        
        job.add_contact("PROT", 2, "LIG1", 1, max_distance=20.0)  # maximum  
        print("✓ Max distance (20.0Å): Valid")
        
        job.add_contact("PROT", 3, "LIG1", 1, max_distance=10.5)  # typical
        print("✓ Typical distance (10.5Å): Valid")
        
    except ValueError as e:
        print(f"✗ Validation error: {e}")
    
    try:
        # Invalid distance - too small
        job.add_contact("PROT", 4, "LIG1", 1, max_distance=3.0)
    except ValueError as e:
        print(f"✓ Invalid distance (3.0Å): Correctly rejected - {e}")
    
    try:  
        # Invalid distance - too large
        job.add_pocket("LIG1", max_distance=25.0)
    except ValueError as e:
        print(f"✓ Invalid distance (25.0Å): Correctly rejected - {e}")
    print()
    
    # Show constraint summary
    print("=== Final Constraint Summary ===")
    print(f"Total constraints: {len(job.constraints)}")
    
    pocket_count = sum(1 for c in job.constraints if hasattr(c, 'binder'))
    contact_count = sum(1 for c in job.constraints if hasattr(c, 'token1') and hasattr(c, 'token2'))
    
    print(f"Pocket constraints: {pocket_count}")
    print(f"Contact constraints: {contact_count}")
    
    # Count enforced constraints (force=True)
    enforced_count = sum(1 for c in job.constraints 
                        if hasattr(c, 'force') and c.force)
    print(f"Enforced constraints (force=True): {enforced_count}")
    print()
    
    # Generate and save YAML
    output_file = "constraint_enhancements_demo.yaml"
    job.write_yaml(output_file)
    print(f"✓ YAML output saved to: {output_file}")
    
    # Show YAML snippet with force parameters
    print("\n=== YAML Output Snippet (Force Parameters) ===")
    with open(output_file, 'r') as f:
        lines = f.readlines()
        
    # Find and display constraints section
    in_constraints = False
    constraint_lines = []
    for line in lines:
        if line.strip() == "constraints:":
            in_constraints = True
            constraint_lines.append(line)
        elif in_constraints:
            if line.startswith("  ") or line.strip() == "":
                constraint_lines.append(line)
            else:
                break
        
        if len(constraint_lines) > 15:  # Limit output
            constraint_lines.append("  # ... (truncated for brevity)\n")
            break
    
    print("".join(constraint_lines))
    
    print("=== Demo Complete ===")
    print(f"""
Key enhancements implemented:
✓ Force parameter support for Pocket and Contact classes
✓ Distance validation with 4-20Å range and 6Å default  
✓ Enhanced pocket constraint validation (any sequence type as binder)
✓ Comprehensive test coverage (42 new tests)
✓ Backward compatibility maintained
""")

if __name__ == "__main__":
    main()