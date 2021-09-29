/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.primerdesign.primer3;


import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

// TODO: make these type safe

/**
 * Java Enum that has constants for all the Primer3 input tags (constructed semi-automatically from
 * the Primer3 manual) to ensure that there is no mis-spelling of tags.
 */
public enum Primer3InputTag {
    // Sequence input tags (specific to a single design)
    SEQUENCE_EXCLUDED_REGION,
    SEQUENCE_INCLUDED_REGION,
    SEQUENCE_PRIMER_REVCOMP,
    SEQUENCE_FORCE_LEFT_END,
    SEQUENCE_INTERNAL_EXCLUDED_REGION,
    SEQUENCE_QUALITY,
    SEQUENCE_FORCE_LEFT_START,
    SEQUENCE_INTERNAL_OLIGO,
    SEQUENCE_START_CODON_POSITION,
    SEQUENCE_FORCE_RIGHT_END,
    SEQUENCE_OVERLAP_JUNCTION_LIST,
    SEQUENCE_TARGET,
    SEQUENCE_FORCE_RIGHT_START,
    SEQUENCE_PRIMER,
    SEQUENCE_TEMPLATE,
    SEQUENCE_ID,
    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST,

    // Global input tags (persist across designs within a primer3 session)
    PRIMER_DNA_CONC,
    PRIMER_MAX_END_GC,
    PRIMER_PAIR_WT_PRODUCT_SIZE_LT,
    PRIMER_DNTP_CONC,
    PRIMER_MAX_END_STABILITY,
    PRIMER_PAIR_WT_PRODUCT_TM_GT,
    PRIMER_EXPLAIN_FLAG,
    PRIMER_MAX_GC,
    PRIMER_PAIR_WT_PRODUCT_TM_LT,
    PRIMER_FIRST_BASE_INDEX,
    PRIMER_MAX_HAIRPIN_TH,
    PRIMER_PAIR_WT_PR_PENALTY,
    PRIMER_GC_CLAMP,
    PRIMER_MAX_LIBRARY_MISPRIMING,
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING,
    PRIMER_INSIDE_PENALTY,
    PRIMER_MAX_NS_ACCEPTED,
    PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH,
    PRIMER_INTERNAL_DNA_CONC,
    PRIMER_MAX_POLY_X,
    PRIMER_PICK_ANYWAY,
    PRIMER_INTERNAL_DNTP_CONC,
    PRIMER_MAX_SELF_ANY,
    PRIMER_PICK_INTERNAL_OLIGO,
    PRIMER_INTERNAL_MAX_GC,
    PRIMER_MAX_SELF_ANY_TH,
    PRIMER_PICK_LEFT_PRIMER,
    PRIMER_INTERNAL_MAX_HAIRPIN_TH,
    PRIMER_MAX_SELF_END,
    PRIMER_PICK_RIGHT_PRIMER,
    PRIMER_INTERNAL_MAX_LIBRARY_MISHYB,
    PRIMER_MAX_SELF_END_TH,
    PRIMER_PRODUCT_MAX_TM,
    PRIMER_INTERNAL_MAX_NS_ACCEPTED,
    PRIMER_MAX_SIZE,
    PRIMER_PRODUCT_MIN_TM,
    PRIMER_INTERNAL_MAX_POLY_X,
    PRIMER_MAX_TEMPLATE_MISPRIMING,
    PRIMER_PRODUCT_OPT_SIZE,
    PRIMER_INTERNAL_MAX_SELF_ANY,
    PRIMER_MAX_TEMPLATE_MISPRIMING_TH,
    PRIMER_PRODUCT_OPT_TM,
    PRIMER_INTERNAL_MAX_SELF_ANY_TH,
    PRIMER_MAX_TM,
    PRIMER_PRODUCT_SIZE_RANGE,
    PRIMER_INTERNAL_MAX_SELF_END,
    PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION,
    PRIMER_QUALITY_RANGE_MAX,
    PRIMER_INTERNAL_MAX_SELF_END_TH,
    PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION,
    PRIMER_QUALITY_RANGE_MIN,
    PRIMER_INTERNAL_MAX_SIZE,
    PRIMER_MIN_END_QUALITY,
    PRIMER_SALT_CORRECTIONS,
    PRIMER_INTERNAL_MAX_TM,
    PRIMER_MIN_GC,
    PRIMER_SALT_DIVALENT,
    PRIMER_INTERNAL_MIN_GC,
    PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE,
    PRIMER_SALT_MONOVALENT,
    PRIMER_INTERNAL_MIN_QUALITY,
    PRIMER_MIN_QUALITY,
    PRIMER_SEQUENCING_ACCURACY,
    PRIMER_INTERNAL_MIN_SIZE,
    PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE,
    PRIMER_SEQUENCING_INTERVAL,
    PRIMER_INTERNAL_MIN_TM,
    PRIMER_MIN_SIZE,
    PRIMER_SEQUENCING_LEAD,
    PRIMER_INTERNAL_MISHYB_LIBRARY,
    PRIMER_MIN_THREE_PRIME_DISTANCE,
    PRIMER_SEQUENCING_SPACING,
    PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME,
    PRIMER_MIN_TM,
    PRIMER_TASK,
    PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME,
    PRIMER_MISPRIMING_LIBRARY,
    PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT,
    PRIMER_INTERNAL_OPT_GC_PERCENT,
    PRIMER_MUST_MATCH_FIVE_PRIME,
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH,
    PRIMER_INTERNAL_OPT_SIZE,
    PRIMER_MUST_MATCH_THREE_PRIME,
    PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT,
    PRIMER_INTERNAL_OPT_TM,
    PRIMER_NUM_RETURN,
    PRIMER_TM_FORMULA,
    PRIMER_INTERNAL_SALT_DIVALENT,
    PRIMER_OPT_GC_PERCENT,
    PRIMER_WT_END_QUAL,
    PRIMER_INTERNAL_SALT_MONOVALENT,
    PRIMER_OPT_SIZE,
    PRIMER_WT_END_STABILITY,
    PRIMER_INTERNAL_WT_END_QUAL,
    PRIMER_OPT_TM,
    PRIMER_WT_GC_PERCENT_GT,
    PRIMER_INTERNAL_WT_GC_PERCENT_GT,
    PRIMER_OUTSIDE_PENALTY,
    PRIMER_WT_GC_PERCENT_LT,
    PRIMER_INTERNAL_WT_GC_PERCENT_LT,
    PRIMER_PAIR_MAX_COMPL_ANY,
    PRIMER_WT_HAIRPIN_TH,
    PRIMER_INTERNAL_WT_HAIRPIN_TH,
    PRIMER_PAIR_MAX_COMPL_ANY_TH,
    PRIMER_WT_LIBRARY_MISPRIMING,
    PRIMER_INTERNAL_WT_LIBRARY_MISHYB,
    PRIMER_PAIR_MAX_COMPL_END,
    PRIMER_WT_NUM_NS,
    PRIMER_INTERNAL_WT_NUM_NS,
    PRIMER_PAIR_MAX_COMPL_END_TH,
    PRIMER_WT_POS_PENALTY,
    PRIMER_INTERNAL_WT_SELF_ANY,
    PRIMER_PAIR_MAX_DIFF_TM,
    PRIMER_WT_SELF_ANY,
    PRIMER_INTERNAL_WT_SELF_ANY_TH,
    PRIMER_PAIR_MAX_LIBRARY_MISPRIMING,
    PRIMER_WT_SELF_ANY_TH,
    PRIMER_INTERNAL_WT_SELF_END,
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING,
    PRIMER_WT_SELF_END,
    PRIMER_INTERNAL_WT_SELF_END_TH,
    PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH,
    PRIMER_WT_SELF_END_TH,
    PRIMER_INTERNAL_WT_SEQ_QUAL,
    PRIMER_PAIR_WT_COMPL_ANY,
    PRIMER_WT_SEQ_QUAL,
    PRIMER_INTERNAL_WT_SIZE_GT,
    PRIMER_PAIR_WT_COMPL_ANY_TH,
    PRIMER_WT_SIZE_GT,
    PRIMER_INTERNAL_WT_SIZE_LT,
    PRIMER_PAIR_WT_COMPL_END,
    PRIMER_WT_SIZE_LT,
    PRIMER_INTERNAL_WT_TM_GT,
    PRIMER_PAIR_WT_COMPL_END_TH,
    PRIMER_WT_TEMPLATE_MISPRIMING,
    PRIMER_INTERNAL_WT_TM_LT,
    PRIMER_PAIR_WT_DIFF_TM,
    PRIMER_WT_TEMPLATE_MISPRIMING_TH,
    PRIMER_LIBERAL_BASE,
    PRIMER_PAIR_WT_IO_PENALTY,
    PRIMER_WT_TM_GT,
    PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS,
    PRIMER_PAIR_WT_LIBRARY_MISPRIMING,
    PRIMER_WT_TM_LT,
    PRIMER_LOWERCASE_MASKING,
    PRIMER_PAIR_WT_PRODUCT_SIZE_GT;

    // The set of all valid input tag names.
    private static final Set<String> names = Arrays.stream(Primer3InputTag.values()).map(x -> x.name()).collect(Collectors.toSet());

    /** Returns true if the provided name represents an input tag, false otherwise. */
    public static boolean includes(final String name) { return names.contains(name); }
}
