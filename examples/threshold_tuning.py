"""
Threshold Validation 

This script validates the default bond displacement threshold by testing multiple 
values across a diverse set of transition state systems.

Test Systems:
- Small molecule reactions (SN2, atropisomer rotation, etc.)
- Organometallic transformations
- Complex multi-bond changes

Metrics:
- Detection Rate: Percentage of expected bonds correctly identified (Recall)
- False Positive Rate: Percentage of detected bonds that are incorrect
- F1 Score: Harmonic mean of precision and recall (balances detection and accuracy)

Usage:
    python examples/tune_thresholds.py
    
Output:
    Console: Progress updates and summary table
    File: threshold_optimization.txt with detailed results
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set, Any
from datetime import datetime

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vib_analysis.api import run_vib_analysis
from vib_analysis import config
from data.expected_results import EXPECTED_RESULTS


def normalize_bond(bond: Tuple[int, int]) -> Tuple[int, int]:
    """Normalize bond tuple to always have lower index first."""
    return tuple(sorted(bond))


def normalize_bonds(bonds: List[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    """Normalize a list of bonds to a set with consistent ordering."""
    return set(normalize_bond(b) for b in bonds)


def test_single_system(
    basename: str,
    bond_threshold: float = 0.4,
) -> Dict[str, Any]:
    """
    Test a single system with specified threshold using the API.
    
    Args:
        basename: System name (without .v000.xyz extension)
        bond_threshold: Bond displacement threshold to test
        
    Returns:
        Dictionary with test results and statistics
    """
    data_dir = Path(__file__).parent / "data"
    xyz_file = data_dir / f"{basename}.v000.xyz"
    
    if not xyz_file.exists():
        raise FileNotFoundError(f"File not found: {xyz_file}")
    
    if basename not in EXPECTED_RESULTS:
        raise ValueError(f"No ground truth for {basename}")
    
    expected = EXPECTED_RESULTS[basename]
    
    # Run analysis using API (includes adaptive threshold logic)
    full_results = run_vib_analysis(
        str(xyz_file),
        bond_threshold=bond_threshold,
        save_trajectory=False,
        print_output=False,
    )
    
    results = full_results['vibrational']
    
    # Compare bonds
    expected_bonds = normalize_bonds(expected.get('bonds', []))
    detected_bonds = normalize_bonds(results['bond_changes'].keys())
    
    bonds_pass = expected_bonds == detected_bonds
    correct_bonds = expected_bonds & detected_bonds
    missing_bonds = expected_bonds - detected_bonds
    extra_bonds = detected_bonds - expected_bonds
    
    # Compare dihedrals (only if expected)
    expected_dihedrals = set(expected.get('dihedrals', []))
    detected_dihedrals = set(results['dihedral_changes'].keys())
    
    dihedrals_pass = expected_dihedrals == detected_dihedrals
    missing_dihedrals = expected_dihedrals - detected_dihedrals
    extra_dihedrals = detected_dihedrals - expected_dihedrals
    
    # Overall pass (only check dihedrals if expected)
    if expected_dihedrals:
        overall_pass = bonds_pass and dihedrals_pass
    else:
        overall_pass = bonds_pass
    
    return {
        'basename': basename,
        'overall_pass': overall_pass,
        'bonds_pass': bonds_pass,
        'expected': expected_bonds,
        'detected': detected_bonds,
        'correct': correct_bonds,
        'missing': missing_bonds,
        'extra': extra_bonds,
        'expected_dihedrals': expected_dihedrals,
        'detected_dihedrals': detected_dihedrals,
        'missing_dihedrals': missing_dihedrals,
        'extra_dihedrals': extra_dihedrals,
        'symbols': results['atom_index_map']
    }


def format_bond(bond: Tuple[int, int], symbols: Dict[int, str]) -> str:
    """Format bond with atom symbols."""
    return f"({bond[0]}, {bond[1]}) [{symbols[bond[0]]}-{symbols[bond[1]]}]"


def write_detailed_results(f, threshold: float, results: List[Dict], stats: Dict):
    """Write detailed per-system results to file."""
    f.write(f"\n{'='*70}\n")
    f.write(f"Testing bond_threshold = {threshold:.2f}\n")
    f.write(f"{'='*70}\n")
    
    for result in results:
        if result['overall_pass']:
            n_bonds = len(result['detected'])
            n_dihedrals = len(result['detected_dihedrals'])
            f.write(f"  ✓ {result['basename']}: Bonds={n_bonds}, Dihedrals={n_dihedrals}\n")
            
            if result['detected']:
                bonds_str = ", ".join(format_bond(b, result['symbols']) 
                                     for b in sorted(result['detected']))
                f.write(f"      Bonds: {bonds_str}\n")
            
            if result['detected_dihedrals']:
                dihedrals_str = ", ".join(str(d) for d in sorted(result['detected_dihedrals']))
                f.write(f"      Dihedrals: {dihedrals_str}\n")
        else:
            basename = result['basename']
            n_exp_bonds = len(result['expected'])
            n_det_bonds = len(result['detected'])
            n_miss_bonds = len(result['missing'])
            n_extra_bonds = len(result['extra'])
            
            n_exp_dih = len(result['expected_dihedrals'])
            n_det_dih = len(result['detected_dihedrals'])
            n_miss_dih = len(result['missing_dihedrals'])
            n_extra_dih = len(result['extra_dihedrals'])
            
            # Build summary line (only show dihedrals if expected)
            summary_parts = []
            if n_exp_bonds > 0 or n_det_bonds > 0:
                summary_parts.append(f"Bonds: Expected={n_exp_bonds}, Detected={n_det_bonds}, Missing={n_miss_bonds}, Extra={n_extra_bonds}")
            if result['expected_dihedrals']:  # Only show if dihedrals are expected
                summary_parts.append(f"Dihedrals: Expected={n_exp_dih}, Detected={n_det_dih}, Missing={n_miss_dih}, Extra={n_extra_dih}")
            
            summary = " | ".join(summary_parts) if summary_parts else "No changes expected or detected"
            f.write(f"  ✗ {basename}: {summary}\n")
            
            # Show bond details
            if result['expected']:
                bonds_str = ", ".join(format_bond(b, result['symbols']) 
                                     for b in sorted(result['expected']))
                f.write(f"      Expected Bonds: {bonds_str}\n")
            if result['detected']:
                bonds_str = ", ".join(format_bond(b, result['symbols']) 
                                     for b in sorted(result['detected']))
                f.write(f"      Detected Bonds: {bonds_str}\n")
            if result['missing']:
                bonds_str = ", ".join(format_bond(b, result['symbols']) 
                                     for b in sorted(result['missing']))
                f.write(f"      Missing Bonds: {bonds_str}\n")
            if result['extra']:
                bonds_str = ", ".join(format_bond(b, result['symbols']) 
                                     for b in sorted(result['extra']))
                f.write(f"      Extra Bonds: {bonds_str}\n")
            
            # Show dihedral details (only if expected)
            if result['expected_dihedrals']:
                if result['expected_dihedrals']:
                    dihedrals_str = ", ".join(str(d) for d in sorted(result['expected_dihedrals']))
                    f.write(f"      Expected Dihedrals: {dihedrals_str}\n")
                if result['detected_dihedrals']:
                    dihedrals_str = ", ".join(str(d) for d in sorted(result['detected_dihedrals']))
                    f.write(f"      Detected Dihedrals: {dihedrals_str}\n")
                if result['missing_dihedrals']:
                    dihedrals_str = ", ".join(str(d) for d in sorted(result['missing_dihedrals']))
                    f.write(f"      Missing Dihedrals: {dihedrals_str}\n")
                if result['extra_dihedrals']:
                    dihedrals_str = ", ".join(str(d) for d in sorted(result['extra_dihedrals']))
                    f.write(f"      Extra Dihedrals: {dihedrals_str}\n")
    
    # Write statistics
    f.write(f"\n  Statistics:\n")
    f.write(f"    Detection Rate: {stats['correct']}/{stats['expected']} bonds ({stats['detection_rate']:.1f}%)\n")
    f.write(f"    False Positives: {stats['extra']} bonds ({stats['false_positive_rate']:.1f}%)\n")
    f.write(f"    F1 Score: {stats['f1_score']:.1f}%\n")
    f.write(f"    Perfect Detection: {stats['passed']}/{stats['total']} systems ({stats['pass_rate']:.1f}%)\n")


def plot_threshold_results(sweep_results: Dict[float, Dict], output_dir: Path):
    """
    Plot threshold validation results showing Detection Rate, False Positive Rate, and F1 Score.
    
    Args:
        sweep_results: Dictionary mapping thresholds to statistics
        output_dir: Directory to save the plot
    """
    # Plotting imports
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Set up global style
    sns.set_theme(style='ticks', context='talk', palette='deep')
    plt.rcParams.update({
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
    })

    plt.rcParams['figure.dpi'] = 100  
    plt.rcParams['savefig.dpi'] = 300 
    # Extract data
    thresholds = sorted(sweep_results.keys())
    detection_rates = [sweep_results[t]['detection_rate'] for t in thresholds]
    false_pos_rates = [sweep_results[t]['false_positive_rate'] for t in thresholds]
    f1_scores = [sweep_results[t]['f1_score'] for t in thresholds]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(5,5))
    
    # Plot lines with markers using custom color ordering
    ax.plot(thresholds, detection_rates, marker='o', linewidth=3, markersize=8, label='Detection', color='teal')
    ax.plot(thresholds, false_pos_rates, marker='s', linewidth=3, markersize=8, label='False\nPositive', color='maroon')
    ax.plot(thresholds, f1_scores, marker='^', linestyle='--', linewidth=3, markersize=8, label='F1 Score', color='darkslateblue')
    
    # Mark current default threshold
    ax.axvline(x=config.BOND_THRESHOLD, color='lightgray', linewidth=20, alpha=0.7, zorder=0)
    
    # Formatting
    ax.set_xlabel('Bond Threshold (Å)')
    ax.set_ylabel('Percentage (%)')
    ax.set_ylim(-5, 110)
    ax.set_xlim(min(thresholds) - 0.02, max(thresholds) + 0.02)

    # Legend
    leg = ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize=14)
    for text in leg.get_texts():
        text.set_ha('center')
    
    ax.set_aspect(1/ax.get_data_ratio())
    
    # Save
    output_path = output_dir / ".." / "images" / "threshold_optimization.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nPlot saved to: {output_path}")


def run_threshold_validation(output_file: str = None):
    """
    Validate bond threshold across multiple values.
    
    Tests thresholds from 0.1 to 0.5 Å and reports optimal value.
    Detailed results written to output_file for transparency.
    """
    # Default output to examples directory
    if output_file is None:
        output_path = Path(__file__).parent / "threshold_optimization.txt"
    else:
        output_path = Path(output_file)
    thresholds = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    
    # Get test systems
    data_dir = Path(__file__).parent / "data"
    xyz_files = sorted(data_dir.glob("*.v000.xyz"))
    basenames = [f.stem.replace('.v000', '') for f in xyz_files 
                 if f.stem.replace('.v000', '') in EXPECTED_RESULTS]
    
    print(f"Validating bond thresholds across {len(basenames)} systems...")
    print(f"Detailed results will be written to: {output_path}\n")
    
    sweep_results = {}
    
    # Open output file with context manager
    with open(output_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("BOND THRESHOLD VALIDATION\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*70 + "\n")
        f.write(f"\nTest systems: {len(basenames)}\n")
        f.write(f"Systems: {', '.join(basenames)}\n")
        
        # Test each threshold
        for idx, threshold in enumerate(thresholds, 1):
            print(f"Testing threshold {idx}/{len(thresholds)}: {threshold:.2f}...", end="\r")
            
            passed = 0
            results = []
            total_expected = 0
            total_correct = 0
            total_extra = 0
            
            for basename in basenames:
                try:
                    result = test_single_system(basename, bond_threshold=threshold)
                    results.append(result)
                    
                    total_expected += len(result['expected'])
                    total_correct += len(result['correct'])
                    total_extra += len(result['extra'])
                    
                    if result['overall_pass']:
                        passed += 1
                        
                except Exception as e:
                    print(f"\nError testing {basename}: {e}")
            
            # Calculate statistics
            total = len(basenames)
            pass_rate = (passed / total * 100) if total > 0 else 0
            detection_rate = (total_correct / total_expected * 100) if total_expected > 0 else 0
            false_pos_rate = (total_extra / (total_correct + total_extra) * 100) if (total_correct + total_extra) > 0 else 0
            
            # Calculate F1 score (harmonic mean of precision and recall)
            precision = (total_correct / (total_correct + total_extra) * 100) if (total_correct + total_extra) > 0 else 0
            recall = detection_rate
            f1_score = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0
            
            stats = {
                'passed': passed,
                'total': total,
                'pass_rate': pass_rate,
                'expected': total_expected,
                'correct': total_correct,
                'extra': total_extra,
                'detection_rate': detection_rate,
                'false_positive_rate': false_pos_rate,
                'f1_score': f1_score
            }
            
            sweep_results[threshold] = stats
            
            # Write detailed results to file
            write_detailed_results(f, threshold, results, stats)
                
        # Write summary to file
        f.write(f"\n{'='*70}\n")
        f.write("SUMMARY\n")
        f.write(f"{'='*70}\n")
        f.write(f"Metric Explanation:\n")
        f.write(f"  Detection Rate: % of expected bonds found (Recall)\n")
        f.write(f"  False Pos Rate: % of detected bonds that are incorrect\n")
        f.write(f"  F1 Score: Harmonic mean of precision and recall\n\n")
        f.write(f"{'Threshold':<12} {'Detection':<15} {'False Pos':<15} {'F1 Score':<12}\n")
        f.write("-" * 70 + "\n")
        
        for threshold in thresholds:
            stats = sweep_results[threshold]
            detection = stats['detection_rate']
            false_pos = stats['false_positive_rate']
            f1 = stats['f1_score']
            
            status = " ← current" if threshold == config.BOND_THRESHOLD else ""
            f.write(f"{threshold:<12.2f} {detection:>5.1f}%          {false_pos:>5.1f}%          {f1:>5.1f}%{status}\n")
       
    # Print summary to console
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print("Metric Explanation:")
    print("  Detection Rate: % of expected bonds found (Recall)")
    print("  False Pos Rate: % of detected bonds that are incorrect")
    print("  F1 Score: Harmonic mean of precision and recall\n")
    print(f"{'Threshold':<12} {'Detection':<15} {'False Pos':<15} {'F1 Score':<12}")
    print("-" * 70)
    
    best_threshold = None
    best_f1 = 0
    
    for threshold in thresholds:
        stats = sweep_results[threshold]
        detection = stats['detection_rate']
        false_pos = stats['false_positive_rate']
        f1 = stats['f1_score']
        
        status = " ← current" if threshold == config.BOND_THRESHOLD else ""
        print(f"{threshold:<12.2f} {detection:>5.1f}%          {false_pos:>5.1f}%          {f1:>5.1f}%{status}")
        
        if f1 > best_f1:
            best_f1 = f1
            best_threshold = threshold
    
    print(f"\n{'='*70}")
    print("VALIDATION RESULT")
    print(f"{'='*70}")
    print(f"Best threshold: {best_threshold:.2f} (F1 Score: {best_f1:.1f}%)")
    
    current_f1 = sweep_results[config.BOND_THRESHOLD]['f1_score']
    if best_threshold == config.BOND_THRESHOLD:
        print(f"✓ Current default ({config.BOND_THRESHOLD}) is optimal!")
    else:
        improvement = best_f1 - current_f1
        print(f"\n!  Recommended: {best_threshold:.2f} (F1: {best_f1:.1f}%)")
        print(f"   Current: {config.BOND_THRESHOLD} (F1: {current_f1:.1f}%)")
        print(f"   Improvement: +{improvement:.1f}%")
    
    print(f"\nDetailed results written to: {output_path}")
    print("="*70)
    
    # Generate plot
    plot_threshold_results(sweep_results, output_path.parent)


if __name__ == "__main__":
    run_threshold_validation()
