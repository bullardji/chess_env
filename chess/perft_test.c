#include "chess.h"
#include <stdio.h>
#include <time.h>
#include <inttypes.h>

// Known perft values for starting position (up to depth 9)
static const uint64_t EXPECTED_PERFT[] = {
    1ULL,                   // 0
    20ULL,                  // 1
    400ULL,                 // 2
    8902ULL,                // 3
    197281ULL,              // 4
    4865609ULL,             // 5
    119060324ULL,           // 6
    3195901860ULL,          // 7
    84998978956ULL,         // 8
    2439530234167ULL        // 9
};

int main(int argc, char** argv) {
    printf("Chess Move Generation Verification (Perft)\n");
    printf("==========================================\n\n");
    
    // Initialize
    init_bitboards();
    printf("✓ Bitboards initialized\n\n");
    
    // Create position
    Position pos;
    pos_set_startpos(&pos);
    
    printf("Testing from starting position:\n");
    printf("rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1\n\n");
    
    printf("Depth   Nodes           Expected        Status\n");
    printf("-----   ------------    ------------    ------\n");
    
    int max_depth = 8; // default
    if (argc > 1) {
        int d = atoi(argv[1]);
        if (d >= 0) max_depth = d;
    }

    // Limit lookups to our EXPECTED_PERFT table size
    int expected_len = (int)(sizeof(EXPECTED_PERFT) / sizeof(EXPECTED_PERFT[0]));

    for (int depth = 0; depth <= max_depth; depth++) {
        clock_t start = clock();
        uint64_t nodes = perft(&pos, depth);
        double elapsed = (double)(clock() - start) / CLOCKS_PER_SEC;
        
        uint64_t expected = (depth < expected_len) ? EXPECTED_PERFT[depth] : nodes;
        const char* status = (depth < expected_len) ? ((nodes == expected) ? "✓ PASS" : "✗ FAIL") : "(no ref)";
        
        printf("  %d     %12" PRIu64 "    %12" PRIu64 "    %s", 
               depth, nodes, expected, status);
        
        if (elapsed > 0.01) {
            printf("  (%.2fs, %.0f knps)", elapsed, nodes / elapsed / 1000.0);
        }
        printf("\n");
        
        if (depth < expected_len && nodes != expected) {
            printf("\nERROR: Perft mismatch at depth %d!\n", depth);
            printf("  Got %" PRIu64 ", expected %" PRIu64 " (diff: %lld)\n", 
                   nodes, expected, (long long)((int64_t)nodes - (int64_t)expected));
            return 1;
        }
    }
    
    printf("\n✓ All perft tests passed! Move generation is correct.\n");
    return 0;
}

