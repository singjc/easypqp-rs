use rustyms::{fragment::FragmentKind, FragmentationModel};
use serde::{Deserialize, Serialize};

pub mod property_prediction;
pub mod tuning_data;
pub mod util;


#[derive(Serialize, Deserialize, Clone)]
pub struct InsilicoPQPSettings {
    /// Precursor charge states to consider
    pub precursor_charge: Vec<u8>,
    /// Maximum charge state of fragment ions to consider. Current MS2 prediction model only supports z=[1,2]
    pub max_fragment_charge: usize,
    /// Minimum number of transitions to consider for a peptide
    pub min_transitions: u8,
    /// Maximum number of transitions to consider for a peptide
    pub max_transitions: u8,
    /// Fragmentation model, (etd/td_etd/ethcd/etcad/eacid/ead/hcd/cid/all/none). See: `[FragmentationModel](https://docs.rs/rustyms/latest/rustyms/model/struct.FragmentationModel.html#method.etd)`
    pub fragmentation_model: String,
    /// Allowed fragment types (default: 'b,y'). Current MS2 prediction model only supports 'b' and 'y'
    pub allowed_fragment_types: Vec<String>,
    /// Scale factor to apply to retention times for output (e.g., 100.0 to convert 0-1 -> 0-100)
    #[serde(default = "InsilicoPQPSettings::default_rt_scale")]
    pub rt_scale: f32,
}

impl InsilicoPQPSettings {
    fn default_rt_scale() -> f32 {
        100.0
    }
}

impl Default for InsilicoPQPSettings {
    fn default() -> Self {
        Self {
            precursor_charge: vec![2, 3],
            max_fragment_charge: 2,
            min_transitions: 6,
            max_transitions: 6,
            fragmentation_model: "cid_hcd".to_string(),
            allowed_fragment_types: vec!["b".to_string(), "y".to_string()],
            rt_scale: 100.0_f32,
        }
    }
    
}

#[derive(Debug, Clone, PartialEq)]
pub struct PeptideProperties {
    pub peptide_index: u32,
    pub retention_time: f32,
    pub precursor: PrecursorProperties,
    pub product: ProductProperties,
}

impl  std::fmt::Display for PeptideProperties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PeptideProperties {{ peptide_index: {}, retention_time: {}, precursor: {:?}, product: {:?} }}", 
            self.peptide_index, self.retention_time, self.precursor, self.product)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrecursorProperties {
    pub peptide_index: u32,
    pub charge: u8,
    pub precursor_mz: f32,
    pub ion_mobility: f64,
}

impl std::fmt::Display for PrecursorProperties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PrecursorProperties {{ peptide_index: {}, charge: {}, precursor_mz: {}, ion_mobility: {:?} }}", 
            self.peptide_index, self.charge, self.precursor_mz, self.ion_mobility)
    }
    
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProductProperties {
    pub peptide_index: u32,
    pub ion_type: Vec<FragmentKind>,
    pub ion_ordinal: Vec<u8>,
    pub charge: Vec<u8>,
    pub product_mz: Vec<f64>,
    pub intensity: Vec<f32>,
}

impl std::fmt::Display for ProductProperties {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ProductProperties {{ peptide_index: {}, ion_type: {:?}, ion_ordinal: {:?}, charge: {:?}, product_mz: {:?}, intensity: {:?} }}", 
            self.peptide_index, self.ion_type, self.ion_ordinal, self.charge, self.product_mz, self.intensity)
    }
    
}

impl ProductProperties {
    /// Returns indices of product ions sorted by intensity (descending)
    pub fn sorted_intensity_indices(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = (0..self.intensity.len()).collect();
        indices.sort_by(|&a, &b| {
            self.intensity[b]
                .partial_cmp(&self.intensity[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        indices
    }

    /// Returns indices of top N most intense product ions
    pub fn top_n_intense_indices(&self, n: usize) -> Vec<usize> {
        let sorted = self.sorted_intensity_indices();
        sorted.into_iter().take(n).collect()
    }
}

/// Selects the fragmentation model based on user input.
/// From https://github.com/rusteomics/mzcore/blob/main/examples/multi-annotator/src/main.rs
pub fn select_model(text: &str, default: &'static FragmentationModel) -> &'static FragmentationModel {
    match text.to_ascii_lowercase().as_str() {
        "etd" => FragmentationModel::etd(),
        "td_etd" => FragmentationModel::td_etd(),
        "ethcd" | "etcad" => FragmentationModel::ethcd(),
        "eacid" => FragmentationModel::eacid(),
        "ead" => FragmentationModel::ead(),
        "hcd" | "cid" => FragmentationModel::cid_hcd(),
        "all" => FragmentationModel::all(),
        "none" => FragmentationModel::none(),
        _ => default,
    }
}

/// Determines if a fragment should be allowed based on current MS2 model capabilities and user preferences.
///
/// # Behavior
/// - Always allows `b` and `y` ions (hard requirement for MS2 model support)
/// - For other fragment types:
///   - Checks against `allowed_types` list
///   - Issues warnings if unsupported types are requested
///   - Still filters them out (only `b`/`y` actually pass)
///
/// # Arguments
/// * `kind` - The fragment type to check
/// * `allowed_types` - User-requested fragment types (e.g., from settings).
///                     Note: Only `b` and `y` are currently supported regardless of this list.
///
/// # Returns
/// `true` if the fragment should be kept (`b` or `y`), `false` otherwise
///
/// # Warning
/// Emits warning messages to stderr when non-`b`/`y` fragments appear in `allowed_types`,
/// as these will still be filtered out due to MS2 model limitations.
///
/// # Example
/// ```ignore
/// // System configuration (would typically come from settings)
/// let allowed_types = vec!["a".into(), "y".into()];
///
/// assert!(is_allowed_fragment(FragmentKind::b, &allowed_types));  // true (always allowed)
/// assert!(is_allowed_fragment(FragmentKind::y, &allowed_types));  // true (always allowed)
/// assert!(!is_allowed_fragment(FragmentKind::a, &allowed_types)); // false (with warning)
/// ```
pub fn is_allowed_fragment(kind: FragmentKind, allowed_types: &[String]) -> bool {
    // Always allow 'b' and 'y' (hard requirement)
    if matches!(kind, FragmentKind::b | FragmentKind::y) {
        return true;
    }

    // Check if this fragment type was requested by the user
    let fragment_str = kind.to_string();
    let was_requested = allowed_types.iter().any(|t| t == &fragment_str);

    // Warn if user asked for unsupported fragments
    if was_requested {
        eprintln!(
            "Warning: MS2 model only supports 'b' and 'y' ions. \
             Discarding requested '{}' fragment.",
            fragment_str
        );
    }

    false // Only b/y actually pass
}


#[cfg(test)]
mod tests {
 
    use super::*;
    use rustyms::fragment::FragmentKind;

    #[test]
    fn test_insilico_default_settings() {
        let settings = InsilicoPQPSettings::default();
        assert_eq!(settings.precursor_charge, vec![2, 3]);
        assert_eq!(settings.max_fragment_charge, 2);
        assert_eq!(settings.min_transitions, 6);
        assert_eq!(settings.max_transitions, 6);
        assert_eq!(settings.fragmentation_model, "cid_hcd");
        assert_eq!(settings.allowed_fragment_types, vec!["b", "y"]);
    }

    #[test]
    fn test_product_properties_sorting() {
        let props = ProductProperties {
            peptide_index: 1,
            ion_type: vec![FragmentKind::b, FragmentKind::y, FragmentKind::b],
            ion_ordinal: vec![1, 2, 3],
            charge: vec![1, 1, 2],
            product_mz: vec![100.0, 200.0, 150.0],
            intensity: vec![10.0, 50.0, 30.0],
        };

        let sorted = props.sorted_intensity_indices();
        assert_eq!(sorted, vec![1, 2, 0]);

        let top2 = props.top_n_intense_indices(2);
        assert_eq!(top2, vec![1, 2]);
    }

    #[test]
    fn test_fragment_filtering_logic() {
        let allowed = vec!["b".to_string(), "y".to_string(), "a".to_string()];
        assert!(is_allowed_fragment(FragmentKind::b, &allowed));
        assert!(is_allowed_fragment(FragmentKind::y, &allowed));
        assert!(!is_allowed_fragment(FragmentKind::a, &allowed)); // warns but returns false
    }

    #[test]
    fn test_select_model() {
        use rustyms::model::FragmentationModel;

        let default = FragmentationModel::cid_hcd();

        assert_eq!(
            select_model("etd", &default),
            FragmentationModel::etd()
        );
        assert_eq!(
            select_model("CID", &default),
            FragmentationModel::cid_hcd()
        );
        assert_eq!(
            select_model("unknown", &default),
            default
        );
    }   

}