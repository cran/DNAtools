## Computes the P_{m/p/Fm/Fp} for a given locus
RPs <- function(p, t, r = 0, R = 0, k = rep(0, 3)) {
  pp <- prob("all", p, R, r, t)  ## Calls c++ functions for computing the relevant probabilities
  ## Mismatch combinations: AABB, AABC, ABCD, 2*AABR_AB, 2*AARB_BA, 2*ARRB_BA, 2*ABCR_ABC,
  ## 2*ABRC_CAB
  mismatch <- pp[c("AABB", "AABC", "ABCD", "AABR_AB", "AABR_AB", "AARB_BA", "AARB_BA", "ARRB_BA", 
    "ARRB_BA", "ABCR_ABC", "ABCR_ABC", "ABRC_CAB", "ABRC_CAB")]
  ## Partial combinations: RARB, AAAB, ABAC, 2*AABR_BA, 2*AARB_AB, BAAR, ABRA, ABBR, BARB,
  ## 2*ARRB_AB ARBR_AB, ARBR_BA, AAAR, 2*ABCR_ACB, 2*ABCR_CAB, 2*ABRC_ABC, 2*ABRC_ACB
  partial <- pp[c("RARB", "AAAB", "ABAC", "AABR_BA", "AABR_BA", "AARB_AB", "AARB_AB", "BAAR", 
    "ABRA", "ABBR", "BARB", "ARRB_AB", "ARRB_AB", "ARBR_AB", "ARBR_BA", "AAAR", "ABCR_ACB", 
    "ABCR_ACB", "ABCR_CAB", "ABCR_CAB", "ABRC_ABC", "ABRC_ABC", "ABRC_ACB", "ABRC_ACB")]
  ## Match combinations: AARR, ABRR, ARAR, ARRA, ARRR, RARA, RARR, RRRR, AAAA, ABAB, ABAR,
  ## BARA, BABR, ABRB
  match <- pp[c("AARR", "ABRR", "ARAR", "ARRA", "ARRR", "RARA", "RARR", "RRRR", "AAAA", "ABAB", 
    "ABAR", "BARA", "BABR", "ABRB")]
  c(sum(mismatch), sum(partial), sum(match))
}
