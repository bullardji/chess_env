#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#ifdef ENABLE_RENDERING
#include "raylib.h"
#endif

typedef uint64_t Bitboard;
typedef uint64_t Key;
typedef uint32_t Square;
typedef uint32_t Move;
typedef uint32_t Piece;
typedef int32_t Value;
typedef int32_t Depth;
typedef uint8_t ChessColor;

#define U64 uint64_t

enum {
    SQ_A1, SQ_B1, SQ_C1, SQ_D1, SQ_E1, SQ_F1, SQ_G1, SQ_H1,
    SQ_A2, SQ_B2, SQ_C2, SQ_D2, SQ_E2, SQ_F2, SQ_G2, SQ_H2,
    SQ_A3, SQ_B3, SQ_C3, SQ_D3, SQ_E3, SQ_F3, SQ_G3, SQ_H3,
    SQ_A4, SQ_B4, SQ_C4, SQ_D4, SQ_E4, SQ_F4, SQ_G4, SQ_H4,
    SQ_A5, SQ_B5, SQ_C5, SQ_D5, SQ_E5, SQ_F5, SQ_G5, SQ_H5,
    SQ_A6, SQ_B6, SQ_C6, SQ_D6, SQ_E6, SQ_F6, SQ_G6, SQ_H6,
    SQ_A7, SQ_B7, SQ_C7, SQ_D7, SQ_E7, SQ_F7, SQ_G7, SQ_H7,
    SQ_A8, SQ_B8, SQ_C8, SQ_D8, SQ_E8, SQ_F8, SQ_G8, SQ_H8,
    SQ_NONE = 64
};

enum { PAWN = 1, KNIGHT, BISHOP, ROOK, QUEEN, KING };

enum {
    NO_PIECE = 0,
    W_PAWN = 1, W_KNIGHT, W_BISHOP, W_ROOK, W_QUEEN, W_KING,
    B_PAWN = 9, B_KNIGHT, B_BISHOP, B_ROOK, B_QUEEN, B_KING
};

enum { CHESS_WHITE = 0, CHESS_BLACK = 1 };

enum {
    NO_CASTLING = 0,
    WHITE_OO = 1, WHITE_OOO = 2,
    BLACK_OO = 4, BLACK_OOO = 8,
    WHITE_CASTLING = 3, BLACK_CASTLING = 12,
    ANY_CASTLING = 15
};


enum { NORMAL, PROMOTION, ENPASSANT, CASTLING };

enum {
    NORTH = 8, EAST = 1, SOUTH = -8, WEST = -1,
    NORTH_EAST = 9, SOUTH_EAST = -7,
    NORTH_WEST = 7, SOUTH_WEST = -9
};

enum { RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8 };
enum { FILE_A, FILE_B, FILE_C, FILE_D, FILE_E, FILE_F, FILE_G, FILE_H };


enum {
    VALUE_ZERO = 0,
    VALUE_DRAW = 0,
    VALUE_MATE = 32000,
    VALUE_INFINITE = 32001,
};


#define PAWN_VALUE 100
#define KNIGHT_VALUE 320
#define BISHOP_VALUE 330
#define ROOK_VALUE 500
#define QUEEN_VALUE 900
#define KING_VALUE 20000

#define MOVE_NONE 0
#define MOVE_NULL 65

#define make_move(from, to) ((Move)((to) | ((from) << 6)))
#define make_promotion(from, to, pt) ((Move)((to) | ((from) << 6) | (PROMOTION << 14) | (((pt) - KNIGHT) << 12)))
#define make_enpassant(from, to) ((Move)((to) | ((from) << 6) | (ENPASSANT << 14)))
#define make_castling(from, to) ((Move)((to) | ((from) << 6) | (CASTLING << 14)))

#define from_sq(m) (((m) >> 6) & 0x3f)
#define to_sq(m) ((m) & 0x3f)
#define type_of_m(m) ((m) >> 14)
#define promotion_type(m) ((((m) >> 12) & 3) + KNIGHT)

#define make_square(f, r) ((Square)(((r) << 3) + (f)))
#define file_of(s) ((s) & 7)
#define rank_of(s) ((s) >> 3)
#define make_piece(c, pt) ((Piece)(((c) << 3) + (pt)))
#define type_of_p(p) ((p) & 7)
#define color_of(p) ((p) >> 3)
#define relative_square(c, s) ((Square)((s) ^ ((c) * 56)))
#define relative_rank(c, r) ((r) ^ ((c) * 7))

#define pieces(pos) ((pos)->byTypeBB[0])
#define pieces_p(pos, p) ((pos)->byTypeBB[p])
#define pieces_c(pos, c) ((pos)->byColorBB[c])
#define pieces_cp(pos, c, p) (pieces_p(pos, p) & pieces_c(pos, c))
#define piece_on(pos, s) ((pos)->board[s])
#define is_empty(pos, s) (piece_on(pos, s) == NO_PIECE)

#define MAX_SERVED_MOVES 64

#define FileABB 0x0101010101010101ULL
#define FileBBB (FileABB << 1)
#define FileCBB (FileABB << 2)
#define FileDBB (FileABB << 3)
#define FileEBB (FileABB << 4)
#define FileFBB (FileABB << 5)
#define FileGBB (FileABB << 6)
#define FileHBB (FileABB << 7)

#define Rank1BB 0xFFULL
#define Rank2BB (Rank1BB << 8)
#define Rank3BB (Rank1BB << 16)
#define Rank4BB (Rank1BB << 24)
#define Rank5BB (Rank1BB << 32)
#define Rank6BB (Rank1BB << 40)
#define Rank7BB (Rank1BB << 48)
#define Rank8BB (Rank1BB << 56)

extern Bitboard SquareBB[65];
extern Bitboard FileBB[8];
extern Bitboard RankBB[8];
extern Bitboard PawnAttacks[2][64];
extern Bitboard KnightAttacks[64];
extern Bitboard KingAttacks[64];
extern Bitboard BetweenBB[64][64];
extern Bitboard LineBB[64][64];

extern Bitboard rook_file_attacks[8][256];
extern Bitboard rook_rank_attacks[8][256];

typedef struct {
    Key psq[16][64];
    Key enpassant[8];
    Key castling[16];
    Key side;
} Zobrist;

extern Zobrist zob;

typedef struct {
    Bitboard byTypeBB[7];    // [0]=all, [1-6]=PAWN,KNIGHT,BISHOP,ROOK,QUEEN,KING
    Bitboard byColorBB[2];
    uint8_t board[64];
    uint8_t pieceCount[16];
    ChessColor sideToMove;
    uint8_t castlingRights;
    uint8_t epSquare;
    uint8_t rule50;
    uint16_t gamePly;
    Key key;
    int16_t materialScore;
    int16_t psqtScore;
} Position;

typedef struct {
    Move move;
    int16_t value;
} ExtMove;

// Transposition table entry
typedef struct {
    Key key;
    Move best_move;
    Value value;
    int8_t depth;
    uint8_t flag;  // 0=exact, 1=lower bound, 2=upper bound
} TTEntry;

#define TT_SIZE 1048576  // 1M entries (~16MB)
#define TT_EXACT 0
#define TT_LOWER 1
#define TT_UPPER 2

static TTEntry transposition_table[TT_SIZE];
static bool tt_initialized = false;

typedef struct {
    ExtMove moves[256];
    int count;
} MoveList;

enum {
    O_BOARD = 0,
    O_SIDE = 768,
    O_CASTLE = 770,
    O_EP = 786,
    O_PICK_PHASE = 851,
    O_SELECTED_PIECE = 853,
    O_VALID_PIECES = 917,
    O_VALID_DESTS = 981,
    OBS_SIZE = 1045
};

typedef struct {
    float score;
    float episode_return;
    float episode_length;
    float wins;
    float losses;
    float draws;
    float timeouts;
    float invalid_action_rate;
    float chess_moves_completed;
    float illegal_moves;
    float avg_legal_moves;
    float n;
} Log;

#ifdef ENABLE_RENDERING
typedef struct {
    Texture2D pieces;
    int cell_size;
} Client;
#else
typedef void Client;
#endif

typedef struct {
    Piece captured;
    uint8_t castlingRights;
    uint8_t epSquare;
    uint8_t rule50;
    int16_t materialScore;
    int16_t psqtScore;
    Key key;
    uint8_t pliesFromNull;
} UndoInfo;

typedef struct {
    Log log;
    Client* client;
    uint8_t* observations;
    int* actions;
    float* rewards;
    unsigned char* terminals;
    
    Position pos;
    MoveList legal_moves;
    int game_result;
    int tick;
    int chess_moves;
    int max_moves;
    int opponent_depth;
    float reward_shaping_weight;
    float reward_draw;
    int render_fps;
    int human_play;
    
    char starting_fen[128];
    
    UndoInfo undo_stack[1024];
    int undo_stack_ptr;
    
    int legal_moves_sum;
    int steps_in_episode;
    int invalid_actions_this_episode;
    int valid_piece_picks;
    
    int pick_phase;
    Square selected_square;
    MoveList valid_destinations;
    float reward_invalid_piece;
    float reward_invalid_move;
    float reward_valid_piece;
    float reward_valid_move;
    
    // Scoreboard tracking
    float ai_score;
    float opponent_score;
    char last_result[32];
} Chess;


static inline Bitboard sq_bb(Square s) {
    return SquareBB[s];
}

static inline int popcount(Bitboard b) {
    return __builtin_popcountll(b);
}

static inline Square lsb(Bitboard b) {
    assert(b);
    return __builtin_ctzll(b);
}

static inline Square pop_lsb(Bitboard* b) {
    Square s = lsb(*b);
    *b &= *b - 1;
    return s;
}

static inline bool more_than_one(Bitboard b) {
    return b & (b - 1);
}

static inline Bitboard shift_bb(int Direction, Bitboard b) {
    return Direction == NORTH ? b << 8
         : Direction == SOUTH ? b >> 8
         : Direction == EAST ? (b & ~FileHBB) << 1
         : Direction == WEST ? (b & ~FileABB) >> 1
         : Direction == NORTH_EAST ? (b & ~FileHBB) << 9
         : Direction == SOUTH_EAST ? (b & ~FileHBB) >> 7
         : Direction == NORTH_WEST ? (b & ~FileABB) << 7
         : Direction == SOUTH_WEST ? (b & ~FileABB) >> 9
         : 0;
}

static inline Bitboard pawn_attacks_bb(ChessColor c, Square s) {
    return PawnAttacks[c][s];
}

static inline Bitboard knight_attacks_bb(Square s) {
    return KnightAttacks[s];
}

static inline Bitboard king_attacks_bb(Square s) {
    return KingAttacks[s];
}


static inline Bitboard rook_attacks_bb(Square s, Bitboard occupied) {
    Bitboard attacks = 0;
    int r = rank_of(s), f = file_of(s);
    for (int rr = r + 1; rr < 8; rr++) {
        Square sq = make_square(f, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int rr = r - 1; rr >= 0; rr--) {
        Square sq = make_square(f, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int ff = f + 1; ff < 8; ff++) {
        Square sq = make_square(ff, r);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int ff = f - 1; ff >= 0; ff--) {
        Square sq = make_square(ff, r);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    return attacks;
}

static inline Bitboard bishop_attacks_bb(Square s, Bitboard occupied) {
    Bitboard attacks = 0;
    int r = rank_of(s), f = file_of(s);
    for (int rr = r + 1, ff = f + 1; rr < 8 && ff < 8; rr++, ff++) {
        Square sq = make_square(ff, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int rr = r - 1, ff = f + 1; rr >= 0 && ff < 8; rr--, ff++) {
        Square sq = make_square(ff, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int rr = r - 1, ff = f - 1; rr >= 0 && ff >= 0; rr--, ff--) {
        Square sq = make_square(ff, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    for (int rr = r + 1, ff = f - 1; rr < 8 && ff >= 0; rr++, ff--) {
        Square sq = make_square(ff, rr);
        attacks |= sq_bb(sq);
        if (occupied & sq_bb(sq)) break;
    }
    return attacks;
}

static inline Bitboard queen_attacks_bb(Square s, Bitboard occupied) {
    return rook_attacks_bb(s, occupied) | bishop_attacks_bb(s, occupied);
}


void init_bitboards(void);
void pos_set(Position* pos, const char* fen);
void pos_set_startpos(Position* pos);

void generate_legal(Position* pos, MoveList* ml, UndoInfo* undo_stack, int* undo_stack_ptr);

void do_move(Position* pos, Move m, UndoInfo* undo_stack, int* undo_stack_ptr);
void undo_move(Position* pos, Move m, UndoInfo* undo_stack, int* undo_stack_ptr);

bool is_check(Position* pos, ChessColor c);
bool is_draw_with_history(Position* pos, UndoInfo* undo_stack, int undo_stack_ptr);

Value evaluate(Position* pos);

uint64_t perft(Position* pos, int depth);

void populate_observations(Chess* env);
void c_reset(Chess* env);
void c_step(Chess* env);
void c_render(Chess* env);
void c_close(Chess* env);

Bitboard SquareBB[65];
Bitboard FileBB[8];
Bitboard RankBB[8];
Bitboard PawnAttacks[2][64];
Bitboard KnightAttacks[64];
Bitboard KingAttacks[64];
Bitboard BetweenBB[64][64];
Bitboard LineBB[64][64];
Zobrist zob;

static bool bitboards_initialized = false;

// Killer moves (2 per ply, max 64 plies)
static Move killer_moves[64][2];

// TT helper functions
static inline void tt_clear(void) {
    memset(transposition_table, 0, sizeof(transposition_table));
    memset(killer_moves, 0, sizeof(killer_moves));
    tt_initialized = true;
}

static inline TTEntry* tt_probe(Key key) {
    return &transposition_table[key % TT_SIZE];
}

static inline void tt_store(Key key, Move move, Value value, int depth, uint8_t flag) {
    TTEntry* entry = tt_probe(key);
    if (entry->depth <= depth) {  // Replace if deeper or equal depth
        entry->key = key;
        entry->best_move = move;
        entry->value = value;
        entry->depth = depth;
        entry->flag = flag;
    }
}

static inline void store_killer(Move move, int ply) {
    if (move != killer_moves[ply][0]) {
        killer_moves[ply][1] = killer_moves[ply][0];
        killer_moves[ply][0] = move;
    }
}

static const int PawnPST[64] = {
      0,   0,   0,   0,   0,   0,   0,   0,
     50,  50,  50,  50,  50,  50,  50,  50,
     10,  10,  20,  30,  30,  20,  10,  10,
      5,   5,  10,  25,  25,  10,   5,   5,
      0,   0,   0,  20,  20,   0,   0,   0,
      5,  -5, -10,   0,   0, -10,  -5,   5,
      5,  10,  10, -20, -20,  10,  10,   5,
      0,   0,   0,   0,   0,   0,   0,   0
};

static const int KnightPST[64] = {
    -50, -40, -30, -30, -30, -30, -40, -50,
    -40, -20,   0,   0,   0,   0, -20, -40,
    -30,   0,  10,  15,  15,  10,   0, -30,
    -30,   5,  15,  20,  20,  15,   5, -30,
    -30,   0,  15,  20,  20,  15,   0, -30,
    -30,   5,  10,  15,  15,  10,   5, -30,
    -40, -20,   0,   5,   5,   0, -20, -40,
    -50, -40, -30, -30, -30, -30, -40, -50
};

static const int BishopPST[64] = {
    -20, -10, -10, -10, -10, -10, -10, -20,
    -10,   0,   0,   0,   0,   0,   0, -10,
    -10,   0,   5,  10,  10,   5,   0, -10,
    -10,   5,   5,  10,  10,   5,   5, -10,
    -10,   0,  10,  10,  10,  10,   0, -10,
    -10,  10,  10,  10,  10,  10,  10, -10,
    -10,   5,   0,   0,   0,   0,   5, -10,
    -20, -10, -10, -10, -10, -10, -10, -20
};

static const int RookPST[64] = {
      0,   0,   0,   0,   0,   0,   0,   0,
      5,  10,  10,  10,  10,  10,  10,   5,
     -5,   0,   0,   0,   0,   0,   0,  -5,
     -5,   0,   0,   0,   0,   0,   0,  -5,
     -5,   0,   0,   0,   0,   0,   0,  -5,
     -5,   0,   0,   0,   0,   0,   0,  -5,
     -5,   0,   0,   0,   0,   0,   0,  -5,
      0,   0,   0,   5,   5,   0,   0,   0
};

static const int QueenPST[64] = {
    -20, -10, -10,  -5,  -5, -10, -10, -20,
    -10,   0,   0,   0,   0,   0,   0, -10,
    -10,   0,   5,   5,   5,   5,   0, -10,
     -5,   0,   5,   5,   5,   5,   0,  -5,
      0,   0,   5,   5,   5,   5,   0,  -5,
    -10,   5,   5,   5,   5,   5,   0, -10,
    -10,   0,   5,   0,   0,   0,   0, -10,
    -20, -10, -10,  -5,  -5, -10, -10, -20
};

static const int KingPST[64] = {
    -30, -40, -40, -50, -50, -40, -40, -30,
    -30, -40, -40, -50, -50, -40, -40, -30,
    -30, -40, -40, -50, -50, -40, -40, -30,
    -30, -40, -40, -50, -50, -40, -40, -30,
    -20, -30, -30, -40, -40, -30, -30, -20,
    -10, -20, -20, -20, -20, -20, -20, -10,
     20,  20,   0,   0,   0,   0,  20,  20,
     20,  30,  10,   0,   0,  10,  30,  20
};


static uint64_t prng_state = 1070372;
static inline uint64_t prng_rand(void) {
    prng_state ^= prng_state >> 12;
    prng_state ^= prng_state << 25;
    prng_state ^= prng_state >> 27;
    return prng_state * 2685821657736338717ULL;
}

void init_bitboards(void) {
    if (bitboards_initialized) return;
    
    for (int c = 0; c < 2; c++) {
        for (int pt = PAWN; pt <= KING; pt++) {
            for (int s = 0; s < 64; s++) {
                zob.psq[make_piece(c, pt)][s] = prng_rand();
            }
        }
    }
    for (int f = 0; f < 8; f++) {
        zob.enpassant[f] = prng_rand();
    }
    for (int cr = 0; cr < 16; cr++) {
        zob.castling[cr] = prng_rand();
    }
    zob.side = prng_rand();
    
    for (int i = 0; i < 64; i++) {
        SquareBB[i] = 1ULL << i;
    }
    SquareBB[64] = 0;
    
    FileBB[0] = FileABB; FileBB[1] = FileBBB; FileBB[2] = FileCBB; FileBB[3] = FileDBB;
    FileBB[4] = FileEBB; FileBB[5] = FileFBB; FileBB[6] = FileGBB; FileBB[7] = FileHBB;
    
    RankBB[0] = Rank1BB; RankBB[1] = Rank2BB; RankBB[2] = Rank3BB; RankBB[3] = Rank4BB;
    RankBB[4] = Rank5BB; RankBB[5] = Rank6BB; RankBB[6] = Rank7BB; RankBB[7] = Rank8BB;
    
    for (int s = 0; s < 64; s++) {
        Bitboard bb = sq_bb(s);
        PawnAttacks[CHESS_WHITE][s] = shift_bb(NORTH_WEST, bb) | shift_bb(NORTH_EAST, bb);
        PawnAttacks[CHESS_BLACK][s] = shift_bb(SOUTH_WEST, bb) | shift_bb(SOUTH_EAST, bb);
    }
    
    int knight_dirs[] = {-17, -15, -10, -6, 6, 10, 15, 17};
    for (int s = 0; s < 64; s++) {
        Bitboard attack = 0;
        int file = file_of(s);
        int rank = rank_of(s);
        
        for (int i = 0; i < 8; i++) {
            int to = s + knight_dirs[i];
            if (to >= 0 && to < 64) {
                int to_file = file_of(to);
                int to_rank = rank_of(to);
                if (abs(to_file - file) <= 2 && abs(to_rank - rank) <= 2) {
                    attack |= sq_bb(to);
                }
            }
        }
        KnightAttacks[s] = attack;
    }
    
    int king_dirs[] = {-9, -8, -7, -1, 1, 7, 8, 9};
    for (int s = 0; s < 64; s++) {
        Bitboard attack = 0;
        int file = file_of(s);
        
        for (int i = 0; i < 8; i++) {
            int to = s + king_dirs[i];
            if (to >= 0 && to < 64) {
                int to_file = file_of(to);
                if (abs(to_file - file) <= 1) {
                    attack |= sq_bb(to);
                }
            }
        }
        KingAttacks[s] = attack;
    }
    
    for (int s1 = 0; s1 < 64; s1++) {
        for (int s2 = 0; s2 < 64; s2++) {
            BetweenBB[s1][s2] = 0;
            LineBB[s1][s2] = 0;
            
            if (s1 == s2) continue;
            
            int f1 = file_of(s1), r1 = rank_of(s1);
            int f2 = file_of(s2), r2 = rank_of(s2);
            int df = f2 - f1, dr = r2 - r1;
            
            if (df == 0 || dr == 0 || abs(df) == abs(dr)) {
                int step_f = df == 0 ? 0 : (df > 0 ? 1 : -1);
                int step_r = dr == 0 ? 0 : (dr > 0 ? 1 : -1);
                
                LineBB[s1][s2] = sq_bb(s1) | sq_bb(s2);
                
                int f = f1 + step_f;
                int r = r1 + step_r;
                
                while (f != f2 || r != r2) {
                    Square sq = make_square(f, r);
                    BetweenBB[s1][s2] |= sq_bb(sq);
                    LineBB[s1][s2] |= sq_bb(sq);
                    f += step_f;
                    r += step_r;
                }
            }
        }
    }
    
    bitboards_initialized = true;
}

void pos_set_startpos(Position* pos) {
    memset(pos, 0, sizeof(Position));
    
    const char* fen = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";
    pos_set(pos, fen);
}

void pos_set(Position* pos, const char* fen) {
    memset(pos, 0, sizeof(Position));
    
    int rank = 7, file = 0;
    const char* ptr = fen;
    
    while (*ptr && *ptr != ' ') {
        char c = *ptr++;
        
        if (c == '/') {
            rank--;
            file = 0;
        } else if (c >= '1' && c <= '8') {
            file += c - '0';
        } else {
            Square sq = make_square(file, rank);
            Piece pc = NO_PIECE;
            int pt = 0, color = 0;
            
            switch (c) {
                case 'P': pc = W_PAWN; pt = PAWN; color = CHESS_WHITE; break;
                case 'N': pc = W_KNIGHT; pt = KNIGHT; color = CHESS_WHITE; break;
                case 'B': pc = W_BISHOP; pt = BISHOP; color = CHESS_WHITE; break;
                case 'R': pc = W_ROOK; pt = ROOK; color = CHESS_WHITE; break;
                case 'Q': pc = W_QUEEN; pt = QUEEN; color = CHESS_WHITE; break;
                case 'K': pc = W_KING; pt = KING; color = CHESS_WHITE; break;
                case 'p': pc = B_PAWN; pt = PAWN; color = CHESS_BLACK; break;
                case 'n': pc = B_KNIGHT; pt = KNIGHT; color = CHESS_BLACK; break;
                case 'b': pc = B_BISHOP; pt = BISHOP; color = CHESS_BLACK; break;
                case 'r': pc = B_ROOK; pt = ROOK; color = CHESS_BLACK; break;
                case 'q': pc = B_QUEEN; pt = QUEEN; color = CHESS_BLACK; break;
                case 'k': pc = B_KING; pt = KING; color = CHESS_BLACK; break;
            }
            
            if (pc != NO_PIECE) {
                pos->board[sq] = pc;
                pos->byTypeBB[pt] |= sq_bb(sq);
                pos->byColorBB[color] |= sq_bb(sq);
                pos->byTypeBB[0] |= sq_bb(sq);
                pos->pieceCount[pc]++;
            }
            file++;
        }
    }
    
    if (*ptr == ' ') ptr++;
    
    pos->sideToMove = (*ptr == 'w') ? CHESS_WHITE : CHESS_BLACK;
    ptr += 2;
    
    pos->castlingRights = NO_CASTLING;
    while (*ptr && *ptr != ' ') {
        if (*ptr == 'K') pos->castlingRights |= WHITE_OO;
        else if (*ptr == 'Q') pos->castlingRights |= WHITE_OOO;
        else if (*ptr == 'k') pos->castlingRights |= BLACK_OO;
        else if (*ptr == 'q') pos->castlingRights |= BLACK_OOO;
        ptr++;
    }
    

    if (*ptr == ' ') ptr++;
    
    pos->epSquare = SQ_NONE;
    if (*ptr != '-') {
        int ep_file = ptr[0] - 'a';
        int ep_rank = ptr[1] - '1';
        pos->epSquare = make_square(ep_file, ep_rank);
    }
    
    static const int piece_value_cp[7] = {0, 100, 320, 330, 500, 900, 0};
    pos->materialScore = 0;
    pos->psqtScore = 0;
    
    for (Square sq = SQ_A1; sq <= SQ_H8; sq++) {
        Piece pc = pos->board[sq];
        if (pc == NO_PIECE) continue;
        
        int pt = type_of_p(pc);
        ChessColor c = color_of(pc);
        int sign = (c == CHESS_WHITE) ? 1 : -1;
        
        pos->materialScore += sign * piece_value_cp[pt];
        
        Square s = (c == CHESS_WHITE) ? sq : (sq ^ 56);
        int pst_val = 0;
        switch (pt) {
            case PAWN: pst_val = PawnPST[s]; break;
            case KNIGHT: pst_val = KnightPST[s]; break;
            case BISHOP: pst_val = BishopPST[s]; break;
            case ROOK: pst_val = RookPST[s]; break;
            case QUEEN: pst_val = QueenPST[s]; break;
            case KING: pst_val = KingPST[s]; break;
        }
        pos->psqtScore += sign * pst_val;
    }
    
    pos->key = 0;
    for (Square sq = SQ_A1; sq <= SQ_H8; sq++) {
        Piece pc = pos->board[sq];
        if (pc != NO_PIECE) {
            pos->key ^= zob.psq[pc][sq];
        }
    }
    if (pos->sideToMove == CHESS_BLACK) {
        pos->key ^= zob.side;
    }
    if (pos->castlingRights) {
        pos->key ^= zob.castling[pos->castlingRights];
    }
    if (pos->epSquare != SQ_NONE) {
        pos->key ^= zob.enpassant[file_of(pos->epSquare)];
    }
}


static void add_move(MoveList* ml, Move m) {
    ml->moves[ml->count].move = m;
    ml->moves[ml->count].value = 0;
    ml->count++;
}

static void generate_pawn_moves(Position* pos, MoveList* ml, ChessColor us) {
    ChessColor them = !us;
    int up = (us == CHESS_WHITE) ? NORTH : SOUTH;
    Bitboard rank7 = (us == CHESS_WHITE) ? Rank7BB : Rank2BB;
    Bitboard rank3 = (us == CHESS_WHITE) ? Rank3BB : Rank6BB;
    
    Bitboard pawns = pieces_cp(pos, us, PAWN);
    Bitboard pawnsOn7 = pawns & rank7;
    Bitboard pawnsNotOn7 = pawns & ~rank7;
    
    Bitboard enemies = pieces_c(pos, them);
    Bitboard empty = ~pieces(pos);
    
    Bitboard b1 = shift_bb(up, pawnsNotOn7) & empty;
    Bitboard b2 = shift_bb(up, b1 & rank3) & empty;
    
    while (b1) {
        Square to = pop_lsb(&b1);
        add_move(ml, make_move(to - up, to));
    }
    
    while (b2) {
        Square to = pop_lsb(&b2);
        add_move(ml, make_move(to - up - up, to));
    }
    
    if (pawnsOn7) {
        Bitboard b3 = shift_bb(up, pawnsOn7) & empty;
        while (b3) {
            Square to = pop_lsb(&b3);
            Square from = to - up;
            add_move(ml, make_promotion(from, to, QUEEN));
            add_move(ml, make_promotion(from, to, ROOK));
            add_move(ml, make_promotion(from, to, BISHOP));
            add_move(ml, make_promotion(from, to, KNIGHT));
        }
    }
    
    Bitboard b4 = shift_bb(up + WEST, pawnsNotOn7) & enemies;
    Bitboard b5 = shift_bb(up + EAST, pawnsNotOn7) & enemies;
    
    while (b4) {
        Square to = pop_lsb(&b4);
        add_move(ml, make_move(to - up - WEST, to));
    }
    
    while (b5) {
        Square to = pop_lsb(&b5);
        add_move(ml, make_move(to - up - EAST, to));
    }
    
    if (pawnsOn7) {
        Bitboard b6 = shift_bb(up + WEST, pawnsOn7) & enemies;
        Bitboard b7 = shift_bb(up + EAST, pawnsOn7) & enemies;
        
        while (b6) {
            Square to = pop_lsb(&b6);
            Square from = to - up - WEST;
            add_move(ml, make_promotion(from, to, QUEEN));
            add_move(ml, make_promotion(from, to, ROOK));
            add_move(ml, make_promotion(from, to, BISHOP));
            add_move(ml, make_promotion(from, to, KNIGHT));
        }
        
        while (b7) {
            Square to = pop_lsb(&b7);
            Square from = to - up - EAST;
            add_move(ml, make_promotion(from, to, QUEEN));
            add_move(ml, make_promotion(from, to, ROOK));
            add_move(ml, make_promotion(from, to, BISHOP));
            add_move(ml, make_promotion(from, to, KNIGHT));
        }
    }
    
    if (pos->epSquare != SQ_NONE) {
        Bitboard ep_pawns = pawnsNotOn7 & pawn_attacks_bb(them, pos->epSquare);
        while (ep_pawns) {
            Square from = pop_lsb(&ep_pawns);
            add_move(ml, make_enpassant(from, pos->epSquare));
        }
    }
}

static void generate_piece_moves(Position* pos, MoveList* ml, int pt, ChessColor us) {
    Bitboard pieces_bb = pieces_cp(pos, us, pt);
    Bitboard target = ~pieces_c(pos, us);
    Bitboard occupied = pieces(pos);
    
    while (pieces_bb) {
        Square from = pop_lsb(&pieces_bb);
        Bitboard attacks = 0;
        
        switch (pt) {
            case KNIGHT:
                attacks = knight_attacks_bb(from);
                break;
            case BISHOP:
                attacks = bishop_attacks_bb(from, occupied);
                break;
            case ROOK:
                attacks = rook_attacks_bb(from, occupied);
                break;
            case QUEEN:
                attacks = queen_attacks_bb(from, occupied);
                break;
            case KING:
                attacks = king_attacks_bb(from);
                break;
        }
        
        attacks &= target;
        
        while (attacks) {
            Square to = pop_lsb(&attacks);
            add_move(ml, make_move(from, to));
        }
    }
}

static bool is_square_attacked(Position* pos, Square sq, ChessColor by_color) {
    Bitboard occupied = pieces(pos);
    
    if (pawn_attacks_bb(!by_color, sq) & pieces_cp(pos, by_color, PAWN))
        return true;
    
    if (knight_attacks_bb(sq) & pieces_cp(pos, by_color, KNIGHT))
        return true;
    
    if (bishop_attacks_bb(sq, occupied) & (pieces_cp(pos, by_color, BISHOP) | pieces_cp(pos, by_color, QUEEN)))
        return true;
    
    if (rook_attacks_bb(sq, occupied) & (pieces_cp(pos, by_color, ROOK) | pieces_cp(pos, by_color, QUEEN)))
        return true;
    
    if (king_attacks_bb(sq) & pieces_cp(pos, by_color, KING))
        return true;
    
    return false;
}

 static void generate_castling(Position* pos, MoveList* ml, ChessColor us) {
    Bitboard occupied = pieces(pos);
    
    if (us == CHESS_WHITE) {
        if (pos->castlingRights & WHITE_OO) {
            if (!(occupied & (sq_bb(SQ_F1) | sq_bb(SQ_G1)))) {
                add_move(ml, make_castling(SQ_E1, SQ_G1));
            }
        }
        if (pos->castlingRights & WHITE_OOO) {
            if (!(occupied & (sq_bb(SQ_D1) | sq_bb(SQ_C1) | sq_bb(SQ_B1)))) {
                add_move(ml, make_castling(SQ_E1, SQ_C1));
            }
        }
    } else {
        if (pos->castlingRights & BLACK_OO) {
            if (!(occupied & (sq_bb(SQ_F8) | sq_bb(SQ_G8)))) {
                add_move(ml, make_castling(SQ_E8, SQ_G8));
            }
        }
        if (pos->castlingRights & BLACK_OOO) {
            if (!(occupied & (sq_bb(SQ_D8) | sq_bb(SQ_C8) | sq_bb(SQ_B8)))) {
                add_move(ml, make_castling(SQ_E8, SQ_C8));
            }
        }
    }
}

static Bitboard attackers_to_sq(Position* pos, Square sq, Bitboard occupied) {
    return (pawn_attacks_bb(CHESS_WHITE, sq) & pieces_cp(pos, CHESS_BLACK, PAWN))
         | (pawn_attacks_bb(CHESS_BLACK, sq) & pieces_cp(pos, CHESS_WHITE, PAWN))
         | (knight_attacks_bb(sq) & pieces_p(pos, KNIGHT))
         | (king_attacks_bb(sq) & pieces_p(pos, KING))
         | (bishop_attacks_bb(sq, occupied) & (pieces_p(pos, BISHOP) | pieces_p(pos, QUEEN)))
         | (rook_attacks_bb(sq, occupied) & (pieces_p(pos, ROOK) | pieces_p(pos, QUEEN)));
}

static inline Bitboard all_pawn_attacks(Bitboard pawns, ChessColor c) {
    if (c == CHESS_WHITE) {
        return ((pawns << 7) & ~FileHBB) | ((pawns << 9) & ~FileABB);
    } else {
        return ((pawns >> 7) & ~FileABB) | ((pawns >> 9) & ~FileHBB);
    }
}

bool is_check(Position* pos, ChessColor c) {
    Bitboard king_bb = pieces_cp(pos, c, KING);
    if (!king_bb) return false;
    Square king_sq = lsb(king_bb);
    return (attackers_to_sq(pos, king_sq, pieces(pos)) & pieces_c(pos, !c)) != 0;
}

// Working baseline evaluation - material + PST only
static Value evaluate_simple(Position* pos) {
    return pos->materialScore + pos->psqtScore;
}

// Enhanced evaluation with careful scaling (based on Stockfish 10)
static Value evaluate_enhanced(Position* pos) {
    // Start with material (most important)
    Value score = pos->materialScore + pos->psqtScore;
    
    ChessColor us = pos->sideToMove;
    ChessColor them = !us;
    Bitboard occupied = pieces(pos);
    
    // Game phase for tapered eval (0=endgame, 24=opening)
    int phase = popcount(pieces_p(pos, KNIGHT)) + popcount(pieces_p(pos, BISHOP))
              + 2 * popcount(pieces_p(pos, ROOK)) + 4 * popcount(pieces_p(pos, QUEEN));
    if (phase > 24) phase = 24;
    
    int mg_bonus = 0, eg_bonus = 0;
    
    // Mobility bonuses (from CFish/Stockfish 10)
    static const int knight_mobility_mg[9] = {-62, -53, -12, -4, 3, 13, 22, 28, 33};
    static const int knight_mobility_eg[9] = {-81, -56, -30, -14, 8, 15, 23, 27, 33};
    static const int bishop_mobility_mg[14] = {-48, -20, 16, 26, 38, 51, 55, 63, 63, 68, 81, 81, 91, 98};
    static const int bishop_mobility_eg[14] = {-59, -23, -3, 13, 24, 42, 54, 57, 65, 73, 78, 86, 88, 97};
    static const int rook_mobility_mg[15] = {-58, -27, -15, -10, -5, -2, 9, 16, 30, 29, 32, 38, 46, 48, 58};
    static const int rook_mobility_eg[15] = {-76, -18, 28, 55, 69, 82, 112, 118, 132, 142, 155, 165, 166, 169, 171};
    static const int queen_mobility_mg[28] = {-39, -21, 3, 3, 14, 22, 28, 41, 43, 48, 56, 60, 60, 66, 67, 70, 71, 73, 79, 88, 88, 99, 102, 102, 106, 109, 113, 116};
    static const int queen_mobility_eg[28] = {-36, -15, 8, 18, 34, 54, 61, 73, 79, 92, 94, 104, 113, 120, 123, 126, 133, 136, 140, 143, 148, 166, 170, 175, 184, 191, 206, 212};
    
    // Mobility area: squares not attacked by enemy pawns, not occupied by our pawns/king
    Bitboard our_pawns = pieces_cp(pos, us, PAWN);
    Bitboard their_pawns = pieces_cp(pos, them, PAWN);
    Bitboard our_pawn_attacks = all_pawn_attacks(our_pawns, us);
    Bitboard their_pawn_attacks = all_pawn_attacks(their_pawns, them);
    Bitboard our_mobility_area = ~(their_pawn_attacks | our_pawns | pieces_cp(pos, us, KING));
    Bitboard their_mobility_area = ~(our_pawn_attacks | their_pawns | pieces_cp(pos, them, KING));
    
    // Knight mobility and outposts
    Bitboard knights = pieces_cp(pos, us, KNIGHT);
    while (knights) {
        Square sq = pop_lsb(&knights);
        int mob = popcount(knight_attacks_bb(sq) & our_mobility_area);
        if (mob < 9) { mg_bonus += knight_mobility_mg[mob]; eg_bonus += knight_mobility_eg[mob]; }
        
        // Outpost: knight on strong square (ranks 4-6, not attacked by enemy pawns)
        int rank = (us == CHESS_WHITE) ? rank_of(sq) : (7 - rank_of(sq));
        if (rank >= 3 && rank <= 5 && !(their_pawn_attacks & sq_bb(sq))) {
            bool supported = (our_pawn_attacks & sq_bb(sq)) != 0;
            mg_bonus += supported ? 36 : 22;
            eg_bonus += supported ? 12 : 6;
        }
    }
    knights = pieces_cp(pos, them, KNIGHT);
    while (knights) {
        Square sq = pop_lsb(&knights);
        int mob = popcount(knight_attacks_bb(sq) & their_mobility_area);
        if (mob < 9) { mg_bonus -= knight_mobility_mg[mob]; eg_bonus -= knight_mobility_eg[mob]; }
        
        int rank = (them == CHESS_WHITE) ? rank_of(sq) : (7 - rank_of(sq));
        if (rank >= 3 && rank <= 5 && !(our_pawn_attacks & sq_bb(sq))) {
            bool supported = (their_pawn_attacks & sq_bb(sq)) != 0;
            mg_bonus -= supported ? 36 : 22;
            eg_bonus -= supported ? 12 : 6;
        }
    }
    
    // Bishop mobility and outposts
    Bitboard bishops = pieces_cp(pos, us, BISHOP);
    while (bishops) {
        Square sq = pop_lsb(&bishops);
        int mob = popcount(bishop_attacks_bb(sq, occupied) & our_mobility_area);
        if (mob < 14) { mg_bonus += bishop_mobility_mg[mob]; eg_bonus += bishop_mobility_eg[mob]; }
        
        // Bishop outpost
        int rank = (us == CHESS_WHITE) ? rank_of(sq) : (7 - rank_of(sq));
        if (rank >= 3 && rank <= 5 && !(their_pawn_attacks & sq_bb(sq))) {
            bool supported = (our_pawn_attacks & sq_bb(sq)) != 0;
            mg_bonus += supported ? 15 : 9;
            eg_bonus += supported ? 5 : 2;
        }
    }
    bishops = pieces_cp(pos, them, BISHOP);
    while (bishops) {
        Square sq = pop_lsb(&bishops);
        int mob = popcount(bishop_attacks_bb(sq, occupied) & their_mobility_area);
        if (mob < 14) { mg_bonus -= bishop_mobility_mg[mob]; eg_bonus -= bishop_mobility_eg[mob]; }
        
        int rank = (them == CHESS_WHITE) ? rank_of(sq) : (7 - rank_of(sq));
        if (rank >= 3 && rank <= 5 && !(our_pawn_attacks & sq_bb(sq))) {
            bool supported = (their_pawn_attacks & sq_bb(sq)) != 0;
            mg_bonus -= supported ? 15 : 9;
            eg_bonus -= supported ? 5 : 2;
        }
    }
    
    // Rook mobility + open files
    Bitboard rooks = pieces_cp(pos, us, ROOK);
    while (rooks) {
        Square sq = pop_lsb(&rooks);
        int file = file_of(sq);
        int mob = popcount(rook_attacks_bb(sq, occupied) & our_mobility_area);
        if (mob < 15) { mg_bonus += rook_mobility_mg[mob]; eg_bonus += rook_mobility_eg[mob]; }
        
        bool our_pawn_on_file = (our_pawns & FileBB[file]) != 0;
        bool their_pawn_on_file = (their_pawns & FileBB[file]) != 0;
        if (!our_pawn_on_file) {
            if (!their_pawn_on_file) {
                mg_bonus += 44; eg_bonus += 20; // Open
            } else {
                mg_bonus += 18; eg_bonus += 7;  // Semi-open
            }
        }
        
        // Rook on pawn (weak pawn)
        if (their_pawn_on_file && !our_pawn_on_file) {
            if (their_pawns & FileBB[file] & RankBB[rank_of(sq)]) {
                mg_bonus += 10; eg_bonus += 28;
            }
        }
    }
    rooks = pieces_cp(pos, them, ROOK);
    while (rooks) {
        Square sq = pop_lsb(&rooks);
        int file = file_of(sq);
        int mob = popcount(rook_attacks_bb(sq, occupied) & their_mobility_area);
        if (mob < 15) { mg_bonus -= rook_mobility_mg[mob]; eg_bonus -= rook_mobility_eg[mob]; }
        
        bool our_pawn_on_file = (our_pawns & FileBB[file]) != 0;
        bool their_pawn_on_file = (their_pawns & FileBB[file]) != 0;
        if (!their_pawn_on_file) {
            if (!our_pawn_on_file) {
                mg_bonus -= 44; eg_bonus -= 20;
            } else {
                mg_bonus -= 18; eg_bonus -= 7;
            }
        }
        
        if (our_pawn_on_file && !their_pawn_on_file) {
            if (our_pawns & FileBB[file] & RankBB[rank_of(sq)]) {
                mg_bonus -= 10; eg_bonus -= 28;
            }
        }
    }
    
    // Queen mobility
    Bitboard queens = pieces_cp(pos, us, QUEEN);
    while (queens) {
        Square sq = pop_lsb(&queens);
        int mob = popcount(queen_attacks_bb(sq, occupied) & our_mobility_area);
        if (mob < 28) { mg_bonus += queen_mobility_mg[mob]; eg_bonus += queen_mobility_eg[mob]; }
    }
    queens = pieces_cp(pos, them, QUEEN);
    while (queens) {
        Square sq = pop_lsb(&queens);
        int mob = popcount(queen_attacks_bb(sq, occupied) & their_mobility_area);
        if (mob < 28) { mg_bonus -= queen_mobility_mg[mob]; eg_bonus -= queen_mobility_eg[mob]; }
    }
    
    // King safety with attack weights
    static const int king_attack_weights[8] = {0, 0, 77, 55, 44, 10};
    Bitboard our_king = pieces_cp(pos, us, KING);
    if (our_king) {
        Square ksq = lsb(our_king);
        Bitboard king_ring = king_attacks_bb(ksq);
        int attack_weight = 0, attack_count = 0;
        
        Bitboard att = pieces_cp(pos, them, KNIGHT);
        while (att) { if (knight_attacks_bb(pop_lsb(&att)) & king_ring) { attack_weight += 77; attack_count++; }}
        att = pieces_cp(pos, them, BISHOP);
        while (att) { if (bishop_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 55; attack_count++; }}
        att = pieces_cp(pos, them, ROOK);
        while (att) { if (rook_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 44; attack_count++; }}
        att = pieces_cp(pos, them, QUEEN);
        while (att) { if (queen_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 10; attack_count++; }}
        
        if (attack_count >= 2) mg_bonus -= attack_weight * attack_count / 2;
        
        // King pawn shield (pawns directly in front of king)
        int king_file = file_of(ksq);
        int shield_pawns = 0;
        for (int f = (king_file > 0 ? king_file - 1 : 0); f <= (king_file < 7 ? king_file + 1 : 7); f++) {
            if (our_pawns & FileBB[f]) shield_pawns++;
        }
        mg_bonus += shield_pawns * 15;
    }
    
    Bitboard their_king = pieces_cp(pos, them, KING);
    if (their_king) {
        Square ksq = lsb(their_king);
        Bitboard king_ring = king_attacks_bb(ksq);
        int attack_weight = 0, attack_count = 0;
        
        Bitboard att = pieces_cp(pos, us, KNIGHT);
        while (att) { if (knight_attacks_bb(pop_lsb(&att)) & king_ring) { attack_weight += 77; attack_count++; }}
        att = pieces_cp(pos, us, BISHOP);
        while (att) { if (bishop_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 55; attack_count++; }}
        att = pieces_cp(pos, us, ROOK);
        while (att) { if (rook_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 44; attack_count++; }}
        att = pieces_cp(pos, us, QUEEN);
        while (att) { if (queen_attacks_bb(pop_lsb(&att), occupied) & king_ring) { attack_weight += 10; attack_count++; }}
        
        if (attack_count >= 2) mg_bonus += attack_weight * attack_count / 2;
        
        int king_file = file_of(ksq);
        int shield_pawns = 0;
        for (int f = (king_file > 0 ? king_file - 1 : 0); f <= (king_file < 7 ? king_file + 1 : 7); f++) {
            if (their_pawns & FileBB[f]) shield_pawns++;
        }
        mg_bonus -= shield_pawns * 15;
    }
    
    // Advanced pawn structure evaluation
    for (int file = 0; file < 8; file++) {
        Bitboard our_pawns_on_file = our_pawns & FileBB[file];
        Bitboard their_pawns_on_file = their_pawns & FileBB[file];
        int our_count = popcount(our_pawns_on_file);
        int their_count = popcount(their_pawns_on_file);
        
        // Doubled pawns
        if (our_count > 1) { mg_bonus -= (our_count - 1) * 11; eg_bonus -= (our_count - 1) * 13; }
        if (their_count > 1) { mg_bonus += (their_count - 1) * 11; eg_bonus += (their_count - 1) * 13; }
        
        // Isolated pawns (no friendly pawns on adjacent files)
        bool our_left = (file > 0) && (our_pawns & FileBB[file - 1]);
        bool our_right = (file < 7) && (our_pawns & FileBB[file + 1]);
        if (our_count > 0 && !our_left && !our_right) {
            mg_bonus -= 11; eg_bonus -= 20; // Isolated pawn penalty
        }
        
        bool their_left = (file > 0) && (their_pawns & FileBB[file - 1]);
        bool their_right = (file < 7) && (their_pawns & FileBB[file + 1]);
        if (their_count > 0 && !their_left && !their_right) {
            mg_bonus += 11; eg_bonus += 20;
        }
    }
    
    // Passed pawns with rank bonuses
    static const int passed_pawn_mg[8] = {0, 5, 12, 10, 57, 163, 271, 0};
    static const int passed_pawn_eg[8] = {0, 18, 23, 31, 62, 167, 250, 0};
    
    Bitboard our_pawn_scan = our_pawns;
    while (our_pawn_scan) {
        Square sq = pop_lsb(&our_pawn_scan);
        int rank = rank_of(sq);
        int file = file_of(sq);
        int our_rank = (us == CHESS_WHITE) ? rank : (7 - rank);
        
        // Check if passed
        Bitboard front_span = 0;
        if (us == CHESS_WHITE) {
            for (int r = rank + 1; r <= 7; r++) {
                front_span |= sq_bb(make_square(file, r));
                if (file > 0) front_span |= sq_bb(make_square(file - 1, r));
                if (file < 7) front_span |= sq_bb(make_square(file + 1, r));
            }
        } else {
            for (int r = rank - 1; r >= 0; r--) {
                front_span |= sq_bb(make_square(file, r));
                if (file > 0) front_span |= sq_bb(make_square(file - 1, r));
                if (file < 7) front_span |= sq_bb(make_square(file + 1, r));
            }
        }
        
        bool is_passed = !(their_pawns & front_span);
        if (is_passed) {
            mg_bonus += passed_pawn_mg[our_rank];
            eg_bonus += passed_pawn_eg[our_rank];
        }
        
        // Pawn phalanx (side-by-side pawns) - stronger than diagonal support
        if ((sq_bb(sq) << 1) & our_pawns) { // Pawn to the left
            mg_bonus += 8 + our_rank * 2;
            eg_bonus += 5 + our_rank;
        }
        
        // Backward pawns: can't advance safely and no friendly pawns behind on adjacent files
        bool can_advance_safely = true;
        Square advance_sq = (us == CHESS_WHITE) ? sq + 8 : sq - 8;
        if (advance_sq >= 0 && advance_sq < 64) {
            if (their_pawn_attacks & sq_bb(advance_sq)) can_advance_safely = false;
        }
        
        bool has_support_behind = false;
        for (int r = (us == CHESS_WHITE ? 0 : 7); r != rank; r += (us == CHESS_WHITE ? 1 : -1)) {
            if (file > 0 && (our_pawns & sq_bb(make_square(file - 1, r)))) has_support_behind = true;
            if (file < 7 && (our_pawns & sq_bb(make_square(file + 1, r)))) has_support_behind = true;
        }
        
        if (!can_advance_safely && !has_support_behind && !is_passed) {
            mg_bonus -= 9; eg_bonus -= 24; // Backward pawn
        }
    }
    
    // Opponent pawns
    Bitboard their_pawn_scan = their_pawns;
    while (their_pawn_scan) {
        Square sq = pop_lsb(&their_pawn_scan);
        int rank = rank_of(sq);
        int file = file_of(sq);
        int their_rank = (them == CHESS_WHITE) ? rank : (7 - rank);
        
        Bitboard front_span = 0;
        if (them == CHESS_WHITE) {
            for (int r = rank + 1; r <= 7; r++) {
                front_span |= sq_bb(make_square(file, r));
                if (file > 0) front_span |= sq_bb(make_square(file - 1, r));
                if (file < 7) front_span |= sq_bb(make_square(file + 1, r));
            }
        } else {
            for (int r = rank - 1; r >= 0; r--) {
                front_span |= sq_bb(make_square(file, r));
                if (file > 0) front_span |= sq_bb(make_square(file - 1, r));
                if (file < 7) front_span |= sq_bb(make_square(file + 1, r));
            }
        }
        
        bool is_passed = !(our_pawns & front_span);
        if (is_passed) {
            mg_bonus -= passed_pawn_mg[their_rank];
            eg_bonus -= passed_pawn_eg[their_rank];
        }
        
        if ((sq_bb(sq) << 1) & their_pawns) {
            mg_bonus -= 8 + their_rank * 2;
            eg_bonus -= 5 + their_rank;
        }
        
        bool can_advance_safely = true;
        Square advance_sq = (them == CHESS_WHITE) ? sq + 8 : sq - 8;
        if (advance_sq >= 0 && advance_sq < 64) {
            if (our_pawn_attacks & sq_bb(advance_sq)) can_advance_safely = false;
        }
        
        bool has_support_behind = false;
        for (int r = (them == CHESS_WHITE ? 0 : 7); r != rank; r += (them == CHESS_WHITE ? 1 : -1)) {
            if (file > 0 && (their_pawns & sq_bb(make_square(file - 1, r)))) has_support_behind = true;
            if (file < 7 && (their_pawns & sq_bb(make_square(file + 1, r)))) has_support_behind = true;
        }
        
        if (!can_advance_safely && !has_support_behind && !is_passed) {
            mg_bonus += 9; eg_bonus += 24;
        }
    }
    
    // Bishop pair bonus
    if (popcount(pieces_cp(pos, us, BISHOP)) >= 2) {
        mg_bonus += 30; eg_bonus += 65;
    }
    if (popcount(pieces_cp(pos, them, BISHOP)) >= 2) {
        mg_bonus -= 30; eg_bonus -= 65;
    }
    
    // Threats: pieces attacking higher value enemy pieces
    Bitboard our_knights = pieces_cp(pos, us, KNIGHT);
    Bitboard our_bishops = pieces_cp(pos, us, BISHOP);
    Bitboard their_queens = pieces_cp(pos, them, QUEEN);
    Bitboard their_rooks = pieces_cp(pos, them, ROOK);
    Bitboard their_bishops = pieces_cp(pos, them, BISHOP);
    Bitboard their_knights = pieces_cp(pos, them, KNIGHT);
    
    // Our minor pieces threatening their pieces
    Bitboard our_minor_attacks = 0;
    Bitboard kn = our_knights;
    while (kn) our_minor_attacks |= knight_attacks_bb(pop_lsb(&kn));
    Bitboard bi = our_bishops;
    while (bi) our_minor_attacks |= bishop_attacks_bb(pop_lsb(&bi), occupied);
    
    if (our_minor_attacks & their_queens) { mg_bonus += 62; eg_bonus += 120; }
    if (our_minor_attacks & their_rooks) { mg_bonus += 68; eg_bonus += 112; }
    if (our_minor_attacks & their_bishops) { mg_bonus += 57; eg_bonus += 44; }
    if (our_minor_attacks & their_knights) { mg_bonus += 39; eg_bonus += 42; }
    
    // Their minor pieces threatening our pieces
    Bitboard their_minor_attacks = 0;
    kn = their_knights;
    while (kn) their_minor_attacks |= knight_attacks_bb(pop_lsb(&kn));
    bi = their_bishops;
    while (bi) their_minor_attacks |= bishop_attacks_bb(pop_lsb(&bi), occupied);
    
    Bitboard our_queens = pieces_cp(pos, us, QUEEN);
    Bitboard our_rooks = pieces_cp(pos, us, ROOK);
    
    if (their_minor_attacks & our_queens) { mg_bonus -= 62; eg_bonus -= 120; }
    if (their_minor_attacks & our_rooks) { mg_bonus -= 68; eg_bonus -= 112; }
    if (their_minor_attacks & our_bishops) { mg_bonus -= 57; eg_bonus -= 44; }
    if (their_minor_attacks & our_knights) { mg_bonus -= 39; eg_bonus -= 42; }
    
    // Rooks threatening pieces
    Bitboard our_rook_attacks = 0;
    Bitboard rk = pieces_cp(pos, us, ROOK);
    while (rk) our_rook_attacks |= rook_attacks_bb(pop_lsb(&rk), occupied);
    
    if (our_rook_attacks & their_queens) { mg_bonus += 51; eg_bonus += 38; }
    if (our_rook_attacks & their_bishops) { mg_bonus += 38; eg_bonus += 71; }
    if (our_rook_attacks & their_knights) { mg_bonus += 38; eg_bonus += 61; }
    
    Bitboard their_rook_attacks = 0;
    rk = pieces_cp(pos, them, ROOK);
    while (rk) their_rook_attacks |= rook_attacks_bb(pop_lsb(&rk), occupied);
    
    if (their_rook_attacks & our_queens) { mg_bonus -= 51; eg_bonus -= 38; }
    if (their_rook_attacks & our_bishops) { mg_bonus -= 38; eg_bonus -= 71; }
    if (their_rook_attacks & our_knights) { mg_bonus -= 38; eg_bonus -= 61; }
    
    // Hanging pieces (undefended pieces)
    Bitboard our_attacks = our_pawn_attacks | our_minor_attacks | our_rook_attacks;
    Bitboard their_pieces = pieces_c(pos, them);
    Bitboard their_defended = their_pawn_attacks;
    
    Bitboard hanging = their_pieces & our_attacks & ~their_defended;
    if (hanging) {
        mg_bonus += popcount(hanging) * 31; // Simplified hanging bonus
        eg_bonus += popcount(hanging) * 17;
    }
    
    Bitboard their_attacks = their_pawn_attacks | their_minor_attacks | their_rook_attacks;
    Bitboard our_pieces = pieces_c(pos, us);
    Bitboard our_defended = our_pawn_attacks;
    
    hanging = our_pieces & their_attacks & ~our_defended;
    if (hanging) {
        mg_bonus -= popcount(hanging) * 31;
        eg_bonus -= popcount(hanging) * 17;
    }
    
    // Tempo bonus (small advantage for side to move)
    mg_bonus += 28; eg_bonus += 28;
    
    // Tapered eval: interpolate between MG and EG based on phase
    score += (mg_bonus * phase + eg_bonus * (24 - phase)) / 24;
    
    return score;
}

// MVV-LVA (Most Valuable Victim - Least Valuable Attacker) move ordering
static int mvv_lva_score(Move m, Position* pos) {
    static const int victim_value[7] = {0, 100, 320, 330, 500, 900, 0};
    static const int attacker_value[7] = {0, 1, 3, 3, 5, 9, 10};
    
    Piece captured = piece_on(pos, to_sq(m));
    Piece attacker = piece_on(pos, from_sq(m));
    
    int score = 0;
    
    // Captures: Victim value * 10 - Attacker value
    if (captured != NO_PIECE) {
        score = victim_value[type_of_p(captured)] * 10 - attacker_value[type_of_p(attacker)];
    }
    
    // Promotions are valuable
    if (type_of_m(m) == PROMOTION) {
        score += victim_value[promotion_type(m)];
    }
    
    // Castling is generally good
    if (type_of_m(m) == CASTLING) {
        score += 50;
    }
    
    return score;
}

// Optimized alpha-beta search with TT, null-move pruning, and LMR
static Value alpha_beta_search(Position* pos, Depth depth, Value alpha, Value beta, int ply, UndoInfo* undo_stack, int* undo_stack_ptr) {
    if (depth <= 0) {
        return evaluate_enhanced(pos);
    }
    
    bool in_check = is_check(pos, pos->sideToMove);
    
    // Transposition table probe
    TTEntry* tt_entry = tt_probe(pos->key);
    if (tt_entry->key == pos->key && tt_entry->depth >= depth) {
        if (tt_entry->flag == TT_EXACT) return tt_entry->value;
        if (tt_entry->flag == TT_LOWER && tt_entry->value >= beta) return tt_entry->value;
        if (tt_entry->flag == TT_UPPER && tt_entry->value <= alpha) return tt_entry->value;
    }
    
    // Null move pruning (if not in check and depth >= 3)
    if (!in_check && depth >= 3 && pos->materialScore + pos->psqtScore > -500) {
        do_move(pos, MOVE_NULL, undo_stack, undo_stack_ptr);
        Value null_score = -alpha_beta_search(pos, depth - 3, -beta, -beta + 1, ply + 1, undo_stack, undo_stack_ptr);
        undo_move(pos, MOVE_NULL, undo_stack, undo_stack_ptr);
        
        if (null_score >= beta) {
            return beta;  // Null move cutoff
        }
    }
    
    MoveList ml;
    generate_legal(pos, &ml, undo_stack, undo_stack_ptr);
    
    if (ml.count == 0) {
        return in_check ? -30000 + ply : 0;  // Checkmate or stalemate
    }
    
    // Enhanced move ordering: TT move > Killers > MVV-LVA
    Move tt_move = (tt_entry->key == pos->key) ? tt_entry->best_move : MOVE_NONE;
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        int score = mvv_lva_score(m, pos);
        
        // Boost TT move
        if (m == tt_move) score += 10000;
        // Boost killer moves
        else if (ply < 64 && (m == killer_moves[ply][0] || m == killer_moves[ply][1])) score += 5000;
        
        ml.moves[i].value = score;
    }
    
    // Sort by score
    for (int i = 0; i < ml.count - 1; i++) {
        for (int j = i + 1; j < ml.count; j++) {
            if (ml.moves[j].value > ml.moves[i].value) {
                ExtMove tmp = ml.moves[i];
                ml.moves[i] = ml.moves[j];
                ml.moves[j] = tmp;
            }
        }
    }
    
    Value best_score = -VALUE_INFINITE;
    Move best_move = ml.moves[0].move;
    uint8_t tt_flag = TT_UPPER;
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        
        // Late Move Reductions (LMR): reduce depth for likely-bad moves
        int new_depth = depth - 1;
        if (i >= 4 && depth >= 3 && !in_check && ml.moves[i].value < 1000) {
            new_depth = depth - 2;  // Search at reduced depth
        }
        
        do_move(pos, m, undo_stack, undo_stack_ptr);
        Value score = -alpha_beta_search(pos, new_depth, -beta, -alpha, ply + 1, undo_stack, undo_stack_ptr);
        
        // Re-search at full depth if LMR fails high
        if (new_depth < depth - 1 && score > alpha) {
            score = -alpha_beta_search(pos, depth - 1, -beta, -alpha, ply + 1, undo_stack, undo_stack_ptr);
        }
        
        undo_move(pos, m, undo_stack, undo_stack_ptr);
        
        if (score > best_score) {
            best_score = score;
            best_move = m;
            
            if (score > alpha) {
                alpha = score;
                tt_flag = TT_EXACT;
                
                if (alpha >= beta) {
                    // Beta cutoff - store killer move
                    if (ply < 64 && piece_on(pos, to_sq(m)) == NO_PIECE) {
                        store_killer(m, ply);
                    }
                    tt_store(pos->key, m, beta, depth, TT_LOWER);
                    return beta;
                }
            }
        }
    }
    
    // Store in TT
    tt_store(pos->key, best_move, best_score, depth, tt_flag);
    return best_score;
}

static Move search_greedy_pseudo_legal(Position* pos, UndoInfo* undo_stack, int* undo_stack_ptr) {
    ChessColor us = pos->sideToMove;
    MoveList ml;
    ml.count = 0;
    generate_pawn_moves(pos, &ml, us);
    generate_piece_moves(pos, &ml, KNIGHT, us);
    generate_piece_moves(pos, &ml, BISHOP, us);
    generate_piece_moves(pos, &ml, ROOK, us);
    generate_piece_moves(pos, &ml, QUEEN, us);
    generate_piece_moves(pos, &ml, KING, us);
    generate_castling(pos, &ml, us);
    
    if (ml.count == 0) return MOVE_NONE;
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        Square to = to_sq(m);
        if (piece_on(pos, to) != NO_PIECE || type_of_m(m) == ENPASSANT) {
            do_move(pos, m, undo_stack, undo_stack_ptr);
            bool legal = !is_check(pos, us);
            undo_move(pos, m, undo_stack, undo_stack_ptr);
            if (legal) return m;
        }
    }
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        Square to = to_sq(m);
        if (piece_on(pos, to) == NO_PIECE && type_of_m(m) != ENPASSANT) {
            do_move(pos, m, undo_stack, undo_stack_ptr);
            bool legal = !is_check(pos, us);
            undo_move(pos, m, undo_stack, undo_stack_ptr);
            if (legal) return m;
        }
    }
    
    return MOVE_NONE;
}

static Move search_random(Position* pos, UndoInfo* undo_stack, int* undo_stack_ptr) {
    MoveList ml;
    generate_legal(pos, &ml, undo_stack, undo_stack_ptr);
    if (ml.count == 0) return MOVE_NONE;
    return ml.moves[rand() % ml.count].move;
}

static Move search_minimax(Position* pos, Depth depth, UndoInfo* undo_stack, int* undo_stack_ptr) {
    // Initialize TT on first use
    if (!tt_initialized) {
        tt_clear();
    }
    
    MoveList ml;
    generate_legal(pos, &ml, undo_stack, undo_stack_ptr);
    
    if (ml.count == 0) return MOVE_NONE;
    
    Move best_move = ml.moves[0].move;
    Value best_score = -32000;
    Value alpha = -32000;
    Value beta = 32000;
    
    // Try TT move first if available
    TTEntry* tt_entry = tt_probe(pos->key);
    Move tt_move = (tt_entry->key == pos->key) ? tt_entry->best_move : MOVE_NONE;
    
    // Enhanced move ordering
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        int score = mvv_lva_score(m, pos);
        if (m == tt_move) score += 10000;  // TT move first
        ml.moves[i].value = score;
    }
    
    // Sort
    for (int i = 0; i < ml.count - 1; i++) {
        for (int j = i + 1; j < ml.count; j++) {
            if (ml.moves[j].value > ml.moves[i].value) {
                ExtMove tmp = ml.moves[i];
                ml.moves[i] = ml.moves[j];
                ml.moves[j] = tmp;
            }
        }
    }
    
    // Evaluate all root moves
    int best_count = 0;
    ExtMove best_moves[256];
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        do_move(pos, m, undo_stack, undo_stack_ptr);
        Value score = -alpha_beta_search(pos, depth - 1, -beta, -alpha, 1, undo_stack, undo_stack_ptr);
        undo_move(pos, m, undo_stack, undo_stack_ptr);
        
        ml.moves[i].value = score;
        
        if (score > best_score) {
            best_score = score;
            best_count = 1;
            best_moves[0] = ml.moves[i];
            if (score > alpha) alpha = score;
        } else if (score == best_score && best_count < 256) {
            // Track all equally good moves
            best_moves[best_count++] = ml.moves[i];
        }
    }
    
    // Add randomness: pick randomly among moves within 50cp of best (wider window)
    int candidate_count = 0;
    ExtMove candidates[256];
    for (int i = 0; i < ml.count; i++) {
        if (ml.moves[i].value >= best_score - 50 && candidate_count < 256) {
            candidates[candidate_count++] = ml.moves[i];
        }
    }
    
    // Pick randomly among best moves
    if (candidate_count > 0) {
        int idx = rand() % candidate_count;
        return candidates[idx].move;
    }
    
    return best_moves[0].move;
}

static Move search_opponent_move(Position* pos, Depth depth, UndoInfo* undo_stack, int* undo_stack_ptr) {
    if (depth < 0) return search_random(pos, undo_stack, undo_stack_ptr);
    if (depth == 0) return search_greedy_pseudo_legal(pos, undo_stack, undo_stack_ptr);
    return search_minimax(pos, depth, undo_stack, undo_stack_ptr);
}

static inline bool is_legal_move(Position* pos, Move m) {
    ChessColor us = pos->sideToMove;
    ChessColor them = (ChessColor)!us;
    int mt = type_of_m(m);
    if (mt == CASTLING) {
        if (is_check(pos, us)) return false;
        Square from = from_sq(m), to = to_sq(m);
        Square mid = (from + to) / 2;
        if (is_square_attacked(pos, mid, them) || is_square_attacked(pos, to, them)) return false;
        return true;
    }
    if (mt == ENPASSANT) {
        Bitboard king_bb = pieces_cp(pos, us, KING);
        if (!king_bb) return false;
        Square ksq = lsb(king_bb);
        Square from = from_sq(m), to = to_sq(m);
        Square capsq = (us == CHESS_WHITE) ? (to - 8) : (to + 8);
        Bitboard occ = pieces(pos) ^ sq_bb(from) ^ sq_bb(capsq) ^ sq_bb(to);
        return (attackers_to_sq(pos, ksq, occ) & pieces_c(pos, them)) == 0;
    }
    UndoInfo u[1]; int p = 0;
    do_move(pos, m, u, &p);
    bool ok = !is_check(pos, us);
    undo_move(pos, m, u, &p);
    return ok;
}

static inline void generate_pseudo_legal(Position* pos, MoveList* ml, ChessColor us) {
    ml->count = 0;
    generate_pawn_moves(pos, ml, us);
    generate_piece_moves(pos, ml, KNIGHT, us);
    generate_piece_moves(pos, ml, BISHOP, us);
    generate_piece_moves(pos, ml, ROOK, us);
    generate_piece_moves(pos, ml, QUEEN, us);
    generate_piece_moves(pos, ml, KING, us);
    generate_castling(pos, ml, us);
}

void generate_legal(Position* pos, MoveList* ml, UndoInfo* undo_stack, int* undo_stack_ptr) {
    MoveList pseudo;
    generate_pseudo_legal(pos, &pseudo, pos->sideToMove);
    ml->count = 0;
    for (int i = 0; i < pseudo.count; i++) {
        if (is_legal_move(pos, pseudo.moves[i].move))
            ml->moves[ml->count++] = pseudo.moves[i];
    }
}


static inline int get_pst_value(Piece pc, Square sq) {
    int pt = type_of_p(pc);
    ChessColor c = color_of(pc);
    Square s = (c == CHESS_WHITE) ? sq : (sq ^ 56);
    
    switch (pt) {
        case PAWN: return PawnPST[s];
        case KNIGHT: return KnightPST[s];
        case BISHOP: return BishopPST[s];
        case ROOK: return RookPST[s];
        case QUEEN: return QueenPST[s];
        case KING: return KingPST[s];
        default: return 0;
    }
}

static inline int get_material_value(Piece pc) {
    static const int piece_value_cp[7] = {0, 100, 320, 330, 500, 900, 0};
    return piece_value_cp[type_of_p(pc)];
}

void do_move(Position* pos, Move m, UndoInfo* undo_stack, int* undo_stack_ptr) {
    // Null move (for null move pruning)
    if (m == MOVE_NULL) {
        undo_stack[*undo_stack_ptr].captured = NO_PIECE;
        undo_stack[*undo_stack_ptr].castlingRights = pos->castlingRights;
        undo_stack[*undo_stack_ptr].epSquare = pos->epSquare;
        undo_stack[*undo_stack_ptr].rule50 = pos->rule50;
        undo_stack[*undo_stack_ptr].materialScore = pos->materialScore;
        undo_stack[*undo_stack_ptr].psqtScore = pos->psqtScore;
        undo_stack[*undo_stack_ptr].key = pos->key;
        undo_stack[*undo_stack_ptr].pliesFromNull = 0;
        (*undo_stack_ptr)++;
        
        if (pos->epSquare != SQ_NONE) {
            pos->key ^= zob.enpassant[file_of(pos->epSquare)];
            pos->epSquare = SQ_NONE;
        }
        pos->sideToMove = !pos->sideToMove;
        pos->key ^= zob.side;
        return;
    }
    
    Square from = from_sq(m);
    Square to = to_sq(m);
    int move_type = type_of_m(m);
    Piece pc = piece_on(pos, from);
    Piece captured = piece_on(pos, to);
    int pt = type_of_p(pc);
    ChessColor us = pos->sideToMove;
    ChessColor them = !us;
    
    undo_stack[*undo_stack_ptr].captured = captured;
    undo_stack[*undo_stack_ptr].castlingRights = pos->castlingRights;
    undo_stack[*undo_stack_ptr].epSquare = pos->epSquare;
    undo_stack[*undo_stack_ptr].rule50 = pos->rule50;
    undo_stack[*undo_stack_ptr].materialScore = pos->materialScore;
    undo_stack[*undo_stack_ptr].psqtScore = pos->psqtScore;
    undo_stack[*undo_stack_ptr].key = pos->key;
    undo_stack[*undo_stack_ptr].pliesFromNull = (*undo_stack_ptr > 0) ? undo_stack[*undo_stack_ptr - 1].pliesFromNull + 1 : 0;
    (*undo_stack_ptr)++;
    
    int sign = (us == CHESS_WHITE) ? 1 : -1;
    int mat_delta = 0, pst_delta = 0;
    
    pst_delta -= sign * get_pst_value(pc, from);
    
    if (captured != NO_PIECE) {
        int cap_sign = (color_of(captured) == CHESS_WHITE) ? 1 : -1;
        mat_delta -= cap_sign * get_material_value(captured);
        pst_delta -= cap_sign * get_pst_value(captured, to);
    }
    
    if (pt == PAWN || captured != NO_PIECE) {
        pos->rule50 = 0;
        undo_stack[*undo_stack_ptr - 1].pliesFromNull = 0;
    }
    else {
        pos->rule50++;
    }
    
    if (pos->epSquare != SQ_NONE) {
        pos->key ^= zob.enpassant[file_of(pos->epSquare)];
    }
    pos->epSquare = SQ_NONE;
    
    if (move_type == CASTLING) {
        pos->key ^= zob.psq[pc][from];
        
        pos->board[from] = NO_PIECE;
        pos->board[to] = pc;
        pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
        pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
        pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);
        pst_delta += sign * get_pst_value(pc, to);
        pos->key ^= zob.psq[pc][to];
        
        Square rook_from, rook_to;
        if (to > from) { 
            rook_from = from + 3;
            rook_to = from + 1;
        } else {
            rook_from = from - 4;
            rook_to = from - 1;
        }
        
        Piece rook = piece_on(pos, rook_from);
        pos->key ^= zob.psq[rook][rook_from];
        pos->board[rook_from] = NO_PIECE;
        pos->board[rook_to] = rook;
        pos->byTypeBB[ROOK] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        pos->byColorBB[us] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        pos->byTypeBB[0] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        pst_delta -= sign * get_pst_value(rook, rook_from);
        pst_delta += sign * get_pst_value(rook, rook_to);
        pos->key ^= zob.psq[rook][rook_to]; 
        
    } else if (move_type == ENPASSANT) {
        pos->key ^= zob.psq[pc][from];
        
        pos->board[from] = NO_PIECE;
        pos->board[to] = pc;
        pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
        pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
        pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);
        pst_delta += sign * get_pst_value(pc, to); 
        pos->key ^= zob.psq[pc][to];
        
        Square cap_sq = to + (us == CHESS_WHITE ? SOUTH : NORTH);
        Piece cap_pawn = piece_on(pos, cap_sq);
        pos->key ^= zob.psq[cap_pawn][cap_sq];
        pos->board[cap_sq] = NO_PIECE;
        pos->byTypeBB[PAWN] ^= sq_bb(cap_sq);
        pos->byColorBB[them] ^= sq_bb(cap_sq);
        pos->byTypeBB[0] ^= sq_bb(cap_sq);
        pos->pieceCount[cap_pawn]--;
        int cap_sign = (color_of(cap_pawn) == CHESS_WHITE) ? 1 : -1;
        mat_delta -= cap_sign * get_material_value(cap_pawn);
        pst_delta -= cap_sign * get_pst_value(cap_pawn, cap_sq);
        
    } else {
        pos->key ^= zob.psq[pc][from];
        
        if (captured != NO_PIECE) {
            pos->key ^= zob.psq[captured][to];
            int cap_pt = type_of_p(captured);
            pos->byTypeBB[cap_pt] ^= sq_bb(to);
            pos->byColorBB[them] ^= sq_bb(to);
            pos->byTypeBB[0] ^= sq_bb(to);
            pos->pieceCount[captured]--;
        }
        
		pos->board[from] = NO_PIECE;
		pos->board[to] = pc;
		pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
		pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
		pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);
        pst_delta += sign * get_pst_value(pc, to);
        pos->key ^= zob.psq[pc][to];
        
        if (move_type == PROMOTION) {
            int promo_pt = promotion_type(m);
            Piece promo_pc = make_piece(us, promo_pt);
            pos->key ^= zob.psq[pc][to];
            pos->board[to] = promo_pc;
            pos->byTypeBB[pt] ^= sq_bb(to);
            pos->byTypeBB[promo_pt] ^= sq_bb(to);
            pos->pieceCount[pc]--;
            pos->pieceCount[promo_pc]++;
            pos->key ^= zob.psq[promo_pc][to];
            mat_delta -= sign * get_material_value(pc);
            mat_delta += sign * get_material_value(promo_pc);
            pst_delta -= sign * get_pst_value(pc, to);
            pst_delta += sign * get_pst_value(promo_pc, to);
        }
        
        if (pt == PAWN) {
            int diff = to - from;
            if (diff == 16 || diff == -16) {
                Square ep_sq = (from + to) / 2;
                if (pawn_attacks_bb(us, ep_sq) & pieces_cp(pos, them, PAWN)) {
                    pos->epSquare = ep_sq;
                    pos->key ^= zob.enpassant[file_of(ep_sq)];
                }
            }
        }
    }
    
    uint8_t old_castling = pos->castlingRights;
    if (pt == KING) {
        pos->castlingRights &= us == CHESS_WHITE ? ~WHITE_CASTLING : ~BLACK_CASTLING;
    }
    if (from == SQ_A1 || to == SQ_A1) pos->castlingRights &= ~WHITE_OOO;
    if (from == SQ_H1 || to == SQ_H1) pos->castlingRights &= ~WHITE_OO;
    if (from == SQ_A8 || to == SQ_A8) pos->castlingRights &= ~BLACK_OOO;
    if (from == SQ_H8 || to == SQ_H8) pos->castlingRights &= ~BLACK_OO;
    
    if (old_castling != pos->castlingRights) {
        pos->key ^= zob.castling[old_castling];
        pos->key ^= zob.castling[pos->castlingRights];
    }
    
    pos->materialScore += mat_delta;
    pos->psqtScore += pst_delta;
    
    pos->sideToMove = them;
    pos->key ^= zob.side;
    pos->gamePly++;
}

void undo_move(Position* pos, Move m, UndoInfo* undo_stack, int* undo_stack_ptr) {
    (*undo_stack_ptr)--;
    UndoInfo* undo = &undo_stack[*undo_stack_ptr];
    
    // Undo null move
    if (m == MOVE_NULL) {
        pos->castlingRights = undo->castlingRights;
        pos->epSquare = undo->epSquare;
        pos->rule50 = undo->rule50;
        pos->materialScore = undo->materialScore;
        pos->psqtScore = undo->psqtScore;
        pos->key = undo->key;
        pos->sideToMove = !pos->sideToMove;
        return;
    }
    
    Square from = from_sq(m);
    Square to = to_sq(m);
    int move_type = type_of_m(m);
    ChessColor us = !pos->sideToMove;
    ChessColor them = pos->sideToMove;
    
    Piece pc = piece_on(pos, to);
    int pt = type_of_p(pc);
    
    pos->castlingRights = undo->castlingRights;
    pos->epSquare = undo->epSquare;
    pos->rule50 = undo->rule50;
    pos->materialScore = undo->materialScore;
    pos->psqtScore = undo->psqtScore;
    pos->key = undo->key;
    pos->sideToMove = us;
    pos->gamePly--;
    
    if (move_type == CASTLING) {
        pos->board[to] = NO_PIECE;
        pos->board[from] = pc;
        pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
        pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
        pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);
        
        Square rook_from, rook_to;
        if (to > from) {
            rook_from = from + 3;
            rook_to = from + 1;
        } else {
            rook_from = from - 4;
            rook_to = from - 1;
        }
        
        Piece rook = piece_on(pos, rook_to);
        pos->board[rook_to] = NO_PIECE;
        pos->board[rook_from] = rook;
        pos->byTypeBB[ROOK] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        pos->byColorBB[us] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        pos->byTypeBB[0] ^= sq_bb(rook_from) ^ sq_bb(rook_to);
        
    } else if (move_type == ENPASSANT) {
        pos->board[to] = NO_PIECE;
        pos->board[from] = pc;
        pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
        pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
        pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);

        
        Square cap_sq = to + (us == CHESS_WHITE ? SOUTH : NORTH);
        Piece cap_pawn = make_piece(them, PAWN);
        pos->board[cap_sq] = cap_pawn;
        pos->byTypeBB[PAWN] ^= sq_bb(cap_sq);
        pos->byColorBB[them] ^= sq_bb(cap_sq);
        pos->byTypeBB[0] ^= sq_bb(cap_sq);
        pos->pieceCount[cap_pawn]++;
        
    } else {
        if (move_type == PROMOTION) {
            int promo_pt = promotion_type(m);
            Piece promo_pc = make_piece(us, promo_pt);
            pc = make_piece(us, PAWN);
            pt = PAWN;
            pos->board[to] = NO_PIECE;
            pos->byTypeBB[promo_pt] ^= sq_bb(to);
            pos->byTypeBB[pt] ^= sq_bb(to);
            pos->pieceCount[promo_pc]--;
            pos->pieceCount[pc]++;
        }
        
        pos->board[to] = undo->captured;
        pos->board[from] = pc;
        pos->byTypeBB[pt] ^= sq_bb(from) ^ sq_bb(to);
        pos->byColorBB[us] ^= sq_bb(from) ^ sq_bb(to);
        
        if (undo->captured != NO_PIECE) {
            int cap_pt = type_of_p(undo->captured);
            pos->byTypeBB[cap_pt] ^= sq_bb(to);
            pos->byColorBB[them] ^= sq_bb(to);
            pos->byTypeBB[0] ^= sq_bb(from);
            pos->pieceCount[undo->captured]++;
        } else {
            pos->byTypeBB[0] ^= sq_bb(from) ^ sq_bb(to);
        }
    }
}

static inline bool is_insufficient_material(const Position* pos) {
    if (pieces_p(pos, PAWN) | pieces_p(pos, ROOK) | pieces_p(pos, QUEEN))
        return false;

    int wN = popcount(pieces_cp(pos, CHESS_WHITE, KNIGHT));
    int bN = popcount(pieces_cp(pos, CHESS_BLACK, KNIGHT));
    int wB = popcount(pieces_cp(pos, CHESS_WHITE, BISHOP));
    int bB = popcount(pieces_cp(pos, CHESS_BLACK, BISHOP));
    int totalMinors = wN + bN + wB + bB;

    if (totalMinors == 0)
        return true;

    if (totalMinors == 1)
        return true;

    if (totalMinors == 2) {
        if ((wN == 2 && wB == 0 && bN == 0 && bB == 0) || (bN == 2 && bB == 0 && wN == 0 && wB == 0))
            return true;
        if ((wN + wB) == 1 && (bN + bB) == 1)
            return true;
    }

    return false;
}

bool is_draw_with_history(Position* pos, UndoInfo* undo_stack, int undo_stack_ptr) {
    if (pos->rule50 >= 100)
        return true;
    
    if (is_insufficient_material(pos))
        return true;
    
    int e = (undo_stack_ptr > 0) ? undo_stack[undo_stack_ptr - 1].pliesFromNull - 4 : -1;
    if (e >= 0) {
        int repetitions = 0;
        for (int i = 4; i <= e + 4; i += 2) {
            if (undo_stack[undo_stack_ptr - 1 - i].key == pos->key) {
                repetitions++;
                if (repetitions >= 2) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

int game_result_with_legal_count(Position* pos, int legal_count, UndoInfo* undo_stack, int undo_stack_ptr) {
    if (legal_count == 0) {
        if (is_check(pos, pos->sideToMove)) {
            return pos->sideToMove == CHESS_WHITE ? 2 : 1;
        } else {
            return 3;
        }
    }
    if (is_draw_with_history(pos, undo_stack, undo_stack_ptr)) {
        return 3;
    }
    return 0;
}


Value evaluate(Position* pos) {
    return evaluate_enhanced(pos);
}


static Value quiesce(Position* pos, Value alpha, Value beta, UndoInfo* undo_stack, int* undo_stack_ptr) {
    Value stand_pat = evaluate_enhanced(pos);
    
    if (stand_pat >= beta)
        return stand_pat;
    
    if (alpha < stand_pat)
        alpha = stand_pat;
    
    MoveList ml;
    ml.count = 0;
    ChessColor us = pos->sideToMove;
    
    generate_pawn_moves(pos, &ml, us);
    generate_piece_moves(pos, &ml, KNIGHT, us);
    generate_piece_moves(pos, &ml, BISHOP, us);
    generate_piece_moves(pos, &ml, ROOK, us);
    generate_piece_moves(pos, &ml, QUEEN, us);
    generate_piece_moves(pos, &ml, KING, us);
    
    UndoInfo local_undo[256];
    int local_ptr = 0;
    
    int capture_count = 0;
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        if (piece_on(pos, to_sq(m)) != NO_PIECE) {
            do_move(pos, m, local_undo, &local_ptr);
            bool legal = !is_check(pos, us);
            undo_move(pos, m, local_undo, &local_ptr);
            
            if (legal) {
                ml.moves[capture_count++] = ml.moves[i];
            }
        }
    }
    ml.count = capture_count;
    
    // Order captures with MVV-LVA
    for (int i = 0; i < ml.count; i++) {
        ml.moves[i].value = mvv_lva_score(ml.moves[i].move, pos);
    }
    for (int i = 0; i < ml.count - 1; i++) {
        for (int j = i + 1; j < ml.count; j++) {
            if (ml.moves[j].value > ml.moves[i].value) {
                ExtMove tmp = ml.moves[i];
                ml.moves[i] = ml.moves[j];
                ml.moves[j] = tmp;
            }
        }
    }
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        do_move(pos, m, local_undo, &local_ptr);
        Value score = -quiesce(pos, -beta, -alpha, local_undo, &local_ptr);
        undo_move(pos, m, local_undo, &local_ptr);
        
        if (score >= beta)
            return beta;
        if (score > alpha)
            alpha = score;
    }
    
    return alpha;
}

static Value alpha_beta(Position* pos, Depth depth, Value alpha, Value beta, UndoInfo* undo_stack, int* undo_stack_ptr) {
    if (depth <= 0)
        return quiesce(pos, alpha, beta, undo_stack, undo_stack_ptr);
    
    MoveList ml;
    generate_legal(pos, &ml, undo_stack, undo_stack_ptr);
    
    if (ml.count == 0) {
        if (is_check(pos, pos->sideToMove))
            return -VALUE_MATE + pos->gamePly;
        return VALUE_DRAW;
    }
    
    // MVV-LVA move ordering
    for (int i = 0; i < ml.count; i++) {
        ml.moves[i].value = mvv_lva_score(ml.moves[i].move, pos);
    }
    
    // Sort moves by score
    for (int i = 0; i < ml.count - 1; i++) {
        for (int j = i + 1; j < ml.count; j++) {
            if (ml.moves[j].value > ml.moves[i].value) {
                ExtMove temp = ml.moves[i];
                ml.moves[i] = ml.moves[j];
                ml.moves[j] = temp;
            }
        }
    }
    
    Value best_score = -VALUE_INFINITE;
    
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        do_move(pos, m, undo_stack, undo_stack_ptr);
        Value score = -alpha_beta(pos, depth - 1, -beta, -alpha, undo_stack, undo_stack_ptr);
        undo_move(pos, m, undo_stack, undo_stack_ptr);
        
        if (score > best_score) {
            best_score = score;
        }
        
        if (score >= beta)
            return beta;
        
        if (score > alpha)
            alpha = score;
    }
    
    return best_score;
}



uint64_t perft(Position* pos, int depth) {
    if (depth == 0)
        return 1ULL;
    
    UndoInfo local_undo[512];
    int local_ptr = 0;
    
    MoveList ml;
    generate_legal(pos, &ml, local_undo, &local_ptr);
    
    uint64_t nodes = 0;
    for (int i = 0; i < ml.count; i++) {
        Move m = ml.moves[i].move;
        do_move(pos, m, local_undo, &local_ptr);
        nodes += perft(pos, depth - 1);
        undo_move(pos, m, local_undo, &local_ptr);
    }
    
    return nodes;
}


void populate_observations(Chess* env) {
    uint8_t* obs = env->observations;
    Position* pos = &env->pos;
    memset(obs, 0, OBS_SIZE);
    
    uint8_t* board_planes = obs + O_BOARD;
    
    for (int sq = 0; sq < 64; sq++) {
        Piece p = pos->board[sq];
        if (p == NO_PIECE) continue;
        
        int plane;
        if (p >= B_PAWN) {
            plane = 6 + (p - B_PAWN);
        } else {
            plane = p - 1;
        }
        board_planes[plane * 64 + sq] = 1;
    }
    
    uint8_t* side_onehot = obs + O_SIDE;
    side_onehot[pos->sideToMove] = 1;
    
    uint8_t* castle_onehot = obs + O_CASTLE;
    castle_onehot[pos->castlingRights] = 1;
    
    uint8_t* ep_onehot = obs + O_EP;
    if (pos->epSquare < 64) {
        ep_onehot[pos->epSquare] = 1;
    } else {
        ep_onehot[64] = 1;
    }
    
    uint8_t* phase_onehot = obs + O_PICK_PHASE;
    phase_onehot[env->pick_phase] = 1;
    
    uint8_t* selected_piece_plane = obs + O_SELECTED_PIECE;
    if (env->pick_phase == 1 && env->selected_square != SQ_NONE) {
        selected_piece_plane[env->selected_square] = 1;
    }
    
    uint8_t* valid_pieces = obs + O_VALID_PIECES;
    uint8_t* valid_dests = obs + O_VALID_DESTS;
    
    if (env->pick_phase == 0) {
        for (int i = 0; i < env->legal_moves.count; i++) {
            Square from = from_sq(env->legal_moves.moves[i].move);
            valid_pieces[from] = 1;
        }
    } else {
        for (int i = 0; i < env->valid_destinations.count; i++) {
            Square to = to_sq(env->valid_destinations.moves[i].move);
            valid_dests[to] = 1;
        }
    }
}

void generate_random_fen(char* fen_out) {
    char board[64];
    memset(board, '.', 64);
    
    // Place kings (not adjacent)
    int wk_sq, bk_sq;
    do {
        wk_sq = rand() % 64;
        bk_sq = rand() % 64;
        int wk_rank = wk_sq / 8, wk_file = wk_sq % 8;
        int bk_rank = bk_sq / 8, bk_file = bk_sq % 8;
        int rank_diff = abs(wk_rank - bk_rank);
        int file_diff = abs(wk_file - bk_file);
        if (wk_sq != bk_sq && (rank_diff > 1 || file_diff > 1)) break;
    } while (1);
    
    board[wk_sq] = 'K';
    board[bk_sq] = 'k';
    
    // Add random pieces (up to 15 white, 15 black)
    const char* white_pieces = "QRRNNBBPP";
    const char* black_pieces = "qrrnnbbpp";
    int num_white = rand() % 16;
    int num_black = rand() % 16;
    
    for (int i = 0; i < num_white; i++) {
        int sq, rank;
        char piece;
        do {
            sq = rand() % 64;
            rank = sq / 8;
            piece = white_pieces[rand() % 9];
        } while (board[sq] != '.' || (piece == 'P' && (rank == 0 || rank == 7)));
        board[sq] = piece;
    }
    
    for (int i = 0; i < num_black; i++) {
        int sq, rank;
        char piece;
        do {
            sq = rand() % 64;
            rank = sq / 8;
            piece = black_pieces[rand() % 9];
        } while (board[sq] != '.' || (piece == 'p' && (rank == 0 || rank == 7)));
        board[sq] = piece;
    }
    
    // Convert to FEN
    char* ptr = fen_out;
    for (int rank = 7; rank >= 0; rank--) {
        int empty = 0;
        for (int file = 0; file < 8; file++) {
            char piece = board[rank * 8 + file];
            if (piece == '.') {
                empty++;
            } else {
                if (empty > 0) {
                    *ptr++ = '0' + empty;
                    empty = 0;
                }
                *ptr++ = piece;
            }
        }
        if (empty > 0) *ptr++ = '0' + empty;
        if (rank > 0) *ptr++ = '/';
    }
    strcpy(ptr, " w - - 0 1");
}

void c_reset(Chess* env) {
    env->tick = 0;
    env->chess_moves = 0;
    env->game_result = 0;
    env->terminals[0] = 0;
    env->rewards[0] = 0.0f;
    env->undo_stack_ptr = 0;
    env->legal_moves_sum = 0;
    env->steps_in_episode = 0;
    env->invalid_actions_this_episode = 0;
    env->valid_piece_picks = 0;
    
    env->pick_phase = 0;
    env->selected_square = SQ_NONE;
    env->valid_destinations.count = 0;
    
    // Keep scoreboard across games (don't reset ai_score/opponent_score/last_result)
    
    if (strcmp(env->starting_fen, "random") == 0) {
        char random_fen[128];
        generate_random_fen(random_fen);
        pos_set(&env->pos, random_fen);
    } else {
        pos_set(&env->pos, env->starting_fen);
    }
    
    generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
    populate_observations(env);
}

void c_step(Chess* env) {
    if (env->human_play) {
        env->rewards[0] = 0.0f;
        env->terminals[0] = 0;
        populate_observations(env);
        return;
    }
    
    if (env->terminals[0]) {
        c_reset(env);
        return;
    }
    
    env->rewards[0] = 0.0f;
    env->terminals[0] = 0;
    
    int action = env->actions[0];
    if (action < 0) action = 0;
    if (action >= 64) action = 63;
    
    Square picked_sq = (Square)action;
    
    if (env->pick_phase == 0) {
        env->tick++;
        Piece pc = piece_on(&env->pos, picked_sq);
        ChessColor us = env->pos.sideToMove;
        
        bool has_legal_moves = false;
        if (pc != NO_PIECE && color_of(pc) == us) {
            env->valid_destinations.count = 0;
            for (int i = 0; i < env->legal_moves.count; i++) {
                Move m = env->legal_moves.moves[i].move;
                if (from_sq(m) == picked_sq) {
                    env->valid_destinations.moves[env->valid_destinations.count++] = env->legal_moves.moves[i];
                    has_legal_moves = true;
                }
            }
        }
        
        if (has_legal_moves) {
            env->valid_piece_picks++;
            env->rewards[0] = env->reward_valid_piece;
            env->selected_square = picked_sq;
            env->pick_phase = 1;
            populate_observations(env);
        } else {
            env->rewards[0] = env->reward_invalid_piece;
            env->invalid_actions_this_episode++;
            populate_observations(env);
        }
        return;
    }
    
    env->tick++;
    
    Move chosen_move = MOVE_NONE;
    for (int i = 0; i < env->valid_destinations.count; i++) {
        if (to_sq(env->valid_destinations.moves[i].move) == picked_sq) {
            chosen_move = env->valid_destinations.moves[i].move;
            break;
        }
    }
    
    if (chosen_move == MOVE_NONE) {
        env->rewards[0] = env->reward_invalid_move;
        env->invalid_actions_this_episode++;
        env->pick_phase = 0;
        env->selected_square = SQ_NONE;
        env->valid_destinations.count = 0;
        populate_observations(env);
        return;
    }
    
    env->chess_moves++;
    env->legal_moves_sum += env->legal_moves.count;
    env->steps_in_episode++;
    
    env->pick_phase = 0;
    env->selected_square = SQ_NONE;
    env->valid_destinations.count = 0;
    
    int16_t material_before = env->pos.materialScore;
    
    do_move(&env->pos, chosen_move, env->undo_stack, &env->undo_stack_ptr);
    
    if (is_check(&env->pos, CHESS_WHITE)) {
        env->log.illegal_moves += 1.0f; 
    }
    
    Move opp_move = search_opponent_move(&env->pos, env->opponent_depth, env->undo_stack, &env->undo_stack_ptr);
    
    if (opp_move == MOVE_NONE) {
        env->rewards[0] = is_check(&env->pos, env->pos.sideToMove) ? 1.0f : env->reward_draw;
        env->terminals[0] = 1;
        if (env->rewards[0] > 0.5f) {
            env->log.wins += 1.0f;
            env->ai_score += 1.0f;
            strcpy(env->last_result, "AI Won!");
        } else {
            env->log.draws += 1.0f;
            env->ai_score += 0.5f;
            env->opponent_score += 0.5f;
            strcpy(env->last_result, "Draw");
        }
        env->log.episode_length += env->tick;
        env->log.invalid_action_rate += (env->tick > 0) ? ((float)env->invalid_actions_this_episode / (float)env->tick) : 0.0f;
        env->log.chess_moves_completed += env->chess_moves;
        env->log.episode_return += env->rewards[0];
        env->log.score += env->rewards[0];
        if (env->steps_in_episode > 0) {
            env->log.avg_legal_moves += (float)env->legal_moves_sum / (float)env->steps_in_episode;
        }
        env->log.n += 1.0f;
        generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
        populate_observations(env);
        return;
    }
    
    if (is_draw_with_history(&env->pos, env->undo_stack, env->undo_stack_ptr)) {
        env->rewards[0] = env->reward_draw;
        env->terminals[0] = 1;
        env->log.draws += 1.0f;
        env->ai_score += 0.5f;
        env->opponent_score += 0.5f;
        strcpy(env->last_result, "Draw");
        env->log.episode_length += env->tick;
        env->log.invalid_action_rate += (env->tick > 0) ? ((float)env->invalid_actions_this_episode / (float)env->tick) : 0.0f;
        env->log.chess_moves_completed += env->chess_moves;
        env->log.episode_return += env->rewards[0];
        env->log.score += env->rewards[0];
        if (env->steps_in_episode > 0) {
            env->log.avg_legal_moves += (float)env->legal_moves_sum / (float)env->steps_in_episode;
        }
        env->log.n += 1.0f;
        generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
        populate_observations(env);
        return;
    }
    
    if (env->chess_moves >= env->max_moves) {
        env->rewards[0] = env->reward_draw;
        env->terminals[0] = 1;
        env->log.timeouts += 1.0f;
        env->ai_score += 0.5f;
        env->opponent_score += 0.5f;
        strcpy(env->last_result, "Timeout (Draw)");
        env->log.episode_length += env->tick;
        env->log.invalid_action_rate += (env->tick > 0) ? ((float)env->invalid_actions_this_episode / (float)env->tick) : 0.0f;
        env->log.chess_moves_completed += env->chess_moves;
        env->log.episode_return += env->rewards[0];
        env->log.score += env->rewards[0];
        if (env->steps_in_episode > 0) {
            env->log.avg_legal_moves += (float)env->legal_moves_sum / (float)env->steps_in_episode;
        }
        env->log.n += 1.0f;
        generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
        populate_observations(env);
        return;
    }
    
    do_move(&env->pos, opp_move, env->undo_stack, &env->undo_stack_ptr);
    
    if (is_check(&env->pos, CHESS_BLACK)) {
        env->log.illegal_moves += 1.0f; 
    }
    
    if (env->undo_stack_ptr > 0 && env->undo_stack[env->undo_stack_ptr - 1].pliesFromNull > 99) {
        env->undo_stack[env->undo_stack_ptr - 1].pliesFromNull = 99;
    }
    
    generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
    env->game_result = game_result_with_legal_count(&env->pos, env->legal_moves.count, env->undo_stack, env->undo_stack_ptr);
    
    if (env->game_result) {
        if (env->game_result == 1) {
            env->rewards[0] = 1.0f;
            env->log.wins += 1.0f;
            env->ai_score += 1.0f;
            strcpy(env->last_result, "AI Won!");
        } else if (env->game_result == 2) {
            env->rewards[0] = -1.0f;
            env->log.losses += 1.0f;
            env->opponent_score += 1.0f;
            strcpy(env->last_result, "Opponent Won");
        } else {
            env->rewards[0] = env->reward_draw;
            env->log.draws += 1.0f;
            env->ai_score += 0.5f;
            env->opponent_score += 0.5f;
            strcpy(env->last_result, "Draw");
        }
        env->terminals[0] = 1;
        env->log.episode_length += env->tick;
        env->log.invalid_action_rate += (env->tick > 0) ? ((float)env->invalid_actions_this_episode / (float)env->tick) : 0.0f;
        env->log.chess_moves_completed += env->chess_moves;
        env->log.episode_return += env->rewards[0];
        env->log.score += env->rewards[0];
        if (env->steps_in_episode > 0) {
            env->log.avg_legal_moves += (float)env->legal_moves_sum / (float)env->steps_in_episode;
        }
        env->log.n += 1.0f;
        populate_observations(env);
        return;
    }
    
    env->rewards[0] = env->reward_valid_move;
    
    if (env->reward_shaping_weight > 0.0f) {
        int16_t material_after = env->pos.materialScore;
        float material_delta = (float)(material_after - material_before);
        float normalized_delta = material_delta / 1800.0f;
        normalized_delta = (normalized_delta > 1.0f) ? 1.0f : (normalized_delta < -1.0f) ? -1.0f : normalized_delta;
        env->rewards[0] += env->reward_shaping_weight * normalized_delta;
    }
    
    populate_observations(env);
}

void c_render(Chess* env) {
#ifdef ENABLE_RENDERING
    const int cell_size = 64;
    const int board_size = 8 * cell_size;
    
    if (env->client == NULL) {
        SetConfigFlags(FLAG_MSAA_4X_HINT);
        InitWindow(board_size, board_size + 80, "PufferLib Chess - AI vs Opponent");
        SetTargetFPS(env->render_fps > 0 ? env->render_fps : 30);
        env->client = (Client*)calloc(1, sizeof(Client));
        env->client->cell_size = cell_size;
        
        // Initialize scoreboard on first render
        env->ai_score = 0.0f;
        env->opponent_score = 0.0f;
        strcpy(env->last_result, "Game starting...");
    }
    
    if (IsKeyDown(KEY_ESCAPE)) {
        CloseWindow();
        exit(0);
    }
    
    static int selected_sq = -1;
    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
        Vector2 mp = GetMousePosition();
        int file = (int)(mp.x) / cell_size;
        int rank = 7 - ((int)(mp.y) / cell_size);
        if (file >= 0 && file < 8 && rank >= 0 && rank < 8) {
            int clicked_sq = (int)make_square(file, rank);
            if (selected_sq == -1) {
                if (env->pos.sideToMove == CHESS_WHITE) {
                    Piece pc = piece_on(&env->pos, (Square)clicked_sq);
                    if (pc != NO_PIECE && color_of(pc) == CHESS_WHITE) {
                        bool has_from = false;
                        for (int i = 0; i < env->legal_moves.count; i++) {
                            if ((int)from_sq(env->legal_moves.moves[i].move) == clicked_sq) { has_from = true; break; }
                        }
                        if (has_from) selected_sq = clicked_sq;
                    }
                }
            } else {
                Move chosen = MOVE_NONE;
                for (int i = 0; i < env->legal_moves.count; i++) {
                    Move m = env->legal_moves.moves[i].move;
                    if ((int)from_sq(m) == selected_sq && (int)to_sq(m) == clicked_sq) { chosen = m; break; }
                }
                if (chosen != MOVE_NONE) {
                    do_move(&env->pos, chosen, env->undo_stack, &env->undo_stack_ptr);
                    env->tick++;
                    Move opp = search_opponent_move(&env->pos, env->opponent_depth, env->undo_stack, &env->undo_stack_ptr);
                    if (opp != MOVE_NONE) {
                        do_move(&env->pos, opp, env->undo_stack, &env->undo_stack_ptr);
                        env->tick++;
                    }
                    generate_legal(&env->pos, &env->legal_moves, env->undo_stack, &env->undo_stack_ptr);
                }
                selected_sq = -1;
            }
        }
    }

    BeginDrawing();
    ClearBackground((Color){40, 40, 40, 255});
    
    for (int rank = 0; rank < 8; rank++) {
        for (int file = 0; file < 8; file++) {
            Color square_color = ((rank + file) % 2 == 0) 
                ? (Color){240, 217, 181, 255}
                : (Color){181, 136, 99, 255};
            
            int draw_x = file * cell_size;
            int draw_y = (7 - rank) * cell_size;
            DrawRectangle(draw_x, draw_y, cell_size, cell_size, square_color);

            if (selected_sq != -1) {
                int sel_f = file_of((Square)selected_sq);
                int sel_r = rank_of((Square)selected_sq);
                if (sel_f == file && sel_r == rank) {
                    DrawRectangleLines(draw_x, draw_y, cell_size, cell_size, (Color){255, 215, 0, 255});
                }
                for (int i = 0; i < env->legal_moves.count; i++) {
                    Move m = env->legal_moves.moves[i].move;
                    if ((int)from_sq(m) == selected_sq) {
                        Square to = to_sq(m);
                        int tf = file_of(to);
                        int tr = rank_of(to);
                        if (tf == file && tr == rank) {
                            DrawRectangleLines(draw_x+2, draw_y+2, cell_size-4, cell_size-4, (Color){0, 200, 0, 255});
                        }
                    }
                }
            }
        }
    }
    
    const char* piece_chars[] = {
        "",      // 0 = empty
        "P", "N", "B", "R", "Q", "K",  // 1..6 white
        "", "",                               // 7,8 unused
        "p", "n", "b", "r", "q", "k"   // 9..14 black
    };
    
    for (Square sq = SQ_A1; sq <= SQ_H8; sq++) {
        Piece pc = piece_on(&env->pos, sq);
        if (pc != NO_PIECE) {
            int file = file_of(sq);
            int rank = rank_of(sq);
            int x = file * cell_size + cell_size / 4;
            int y = (7 - rank) * cell_size + cell_size / 8;
            
            Color pc_color = color_of(pc) == CHESS_WHITE 
                ? (Color){255, 255, 255, 255}
                : (Color){0, 0, 0, 255};
            
            DrawText(piece_chars[pc], x, y, cell_size / 2, pc_color);
        }
    }
    
    // Draw scoreboard at bottom
    const int scoreboard_y = board_size + 10;
    char score_text[128];
    snprintf(score_text, sizeof(score_text), "AI: %.1f  Opponent: %.1f", 
             env->ai_score, env->opponent_score);
    DrawText(score_text, 10, scoreboard_y, 20, WHITE);
    
    // Draw last result
    if (env->last_result[0] != '\0') {
        Color result_color = GREEN;
        if (strstr(env->last_result, "Opponent")) result_color = RED;
        else if (strstr(env->last_result, "Draw")) result_color = YELLOW;
        
        DrawText(env->last_result, 10, scoreboard_y + 25, 18, result_color);
    }
    
    // Draw current move count
    char move_text[64];
    snprintf(move_text, sizeof(move_text), "Move: %d", env->chess_moves);
    DrawText(move_text, board_size - 100, scoreboard_y, 18, LIGHTGRAY);
    
    EndDrawing();
#else
    (void)env;  // Unused when rendering disabled
#endif
}

void c_close(Chess* env) {
#ifdef ENABLE_RENDERING
    if (env->client != NULL) {
        if (IsWindowReady()) {
            CloseWindow();
        }
        free(env->client);
        env->client = NULL;
    }
#else
    (void)env;  // Unused when rendering disabled
#endif
}

