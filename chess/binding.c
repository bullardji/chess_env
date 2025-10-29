#include "chess.h"

#define Env Chess
#define MY_SHARED
#include "../env_binding.h"

static PyObject* my_shared(PyObject* self, PyObject* args, PyObject* kwargs) {
    init_bitboards();
    Py_RETURN_NONE;
}

static int my_init(Env *env, PyObject *args, PyObject *kwargs) {
    static int bitboards_initialized = 0;
    if (!bitboards_initialized) {
        init_bitboards();
        bitboards_initialized = 1;
    }
    
    env->max_moves = 500;
    env->opponent_depth = 1;
    env->reward_shaping_weight = 1.0f;
    env->reward_draw = 0.0f;
    env->reward_invalid_piece = -0.01f;
    env->reward_invalid_move = -0.01f;
    env->reward_valid_piece = 0.01f;
    env->reward_valid_move = 0.01f;
    env->client = NULL;
    env->render_fps = 30;
    env->human_play = 0;
    strcpy(env->starting_fen, "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1");
    
    if (kwargs != NULL) {
        PyObject* max_moves_obj = PyDict_GetItemString(kwargs, "max_moves");
        if (max_moves_obj != NULL && PyLong_Check(max_moves_obj)) {
            env->max_moves = (int)PyLong_AsLong(max_moves_obj);
        }
        
        PyObject* opponent_depth_obj = PyDict_GetItemString(kwargs, "opponent_depth");
        if (opponent_depth_obj != NULL && PyLong_Check(opponent_depth_obj)) {
            env->opponent_depth = (int)PyLong_AsLong(opponent_depth_obj);
        }
        
        PyObject* shaping_obj = PyDict_GetItemString(kwargs, "reward_shaping_weight");
        if (shaping_obj != NULL && PyFloat_Check(shaping_obj)) {
            env->reward_shaping_weight = (float)PyFloat_AsDouble(shaping_obj);
        } else if (shaping_obj != NULL && PyLong_Check(shaping_obj)) {
            env->reward_shaping_weight = (float)PyLong_AsDouble(shaping_obj);
        }

        PyObject* reward_draw_obj = PyDict_GetItemString(kwargs, "reward_draw");
        if (reward_draw_obj != NULL && PyFloat_Check(reward_draw_obj)) {
            env->reward_draw = (float)PyFloat_AsDouble(reward_draw_obj);
        } else if (reward_draw_obj != NULL && PyLong_Check(reward_draw_obj)) {
            env->reward_draw = (float)PyLong_AsDouble(reward_draw_obj);
        }

        PyObject* reward_invalid_piece_obj = PyDict_GetItemString(kwargs, "reward_invalid_piece");
        if (reward_invalid_piece_obj != NULL && PyFloat_Check(reward_invalid_piece_obj)) {
            env->reward_invalid_piece = (float)PyFloat_AsDouble(reward_invalid_piece_obj);
        } else if (reward_invalid_piece_obj != NULL && PyLong_Check(reward_invalid_piece_obj)) {
            env->reward_invalid_piece = (float)PyLong_AsDouble(reward_invalid_piece_obj);
        }

        PyObject* reward_invalid_move_obj = PyDict_GetItemString(kwargs, "reward_invalid_move");
        if (reward_invalid_move_obj != NULL && PyFloat_Check(reward_invalid_move_obj)) {
            env->reward_invalid_move = (float)PyFloat_AsDouble(reward_invalid_move_obj);
        } else if (reward_invalid_move_obj != NULL && PyLong_Check(reward_invalid_move_obj)) {
            env->reward_invalid_move = (float)PyLong_AsDouble(reward_invalid_move_obj);
        }

        PyObject* reward_valid_piece_obj = PyDict_GetItemString(kwargs, "reward_valid_piece");
        if (reward_valid_piece_obj != NULL && PyFloat_Check(reward_valid_piece_obj)) {
            env->reward_valid_piece = (float)PyFloat_AsDouble(reward_valid_piece_obj);
        } else if (reward_valid_piece_obj != NULL && PyLong_Check(reward_valid_piece_obj)) {
            env->reward_valid_piece = (float)PyLong_AsDouble(reward_valid_piece_obj);
        }

        PyObject* reward_valid_move_obj = PyDict_GetItemString(kwargs, "reward_valid_move");
        if (reward_valid_move_obj != NULL && PyFloat_Check(reward_valid_move_obj)) {
            env->reward_valid_move = (float)PyFloat_AsDouble(reward_valid_move_obj);
        } else if (reward_valid_move_obj != NULL && PyLong_Check(reward_valid_move_obj)) {
            env->reward_valid_move = (float)PyLong_AsDouble(reward_valid_move_obj);
        }

        PyObject* fps_obj = PyDict_GetItemString(kwargs, "render_fps");
        if (fps_obj != NULL && PyLong_Check(fps_obj)) {
            env->render_fps = (int)PyLong_AsLong(fps_obj);
        }

        PyObject* human_obj = PyDict_GetItemString(kwargs, "human_play");
        if (human_obj != NULL && PyLong_Check(human_obj)) {
            env->human_play = (int)PyLong_AsLong(human_obj);
        }

        PyObject* fen_obj = PyDict_GetItemString(kwargs, "starting_fen");
        if (fen_obj != NULL && PyUnicode_Check(fen_obj)) {
            const char* fen_str = PyUnicode_AsUTF8(fen_obj);
            if (fen_str != NULL) {
                strncpy(env->starting_fen, fen_str, sizeof(env->starting_fen) - 1);
                env->starting_fen[sizeof(env->starting_fen) - 1] = '\0';
            }
        }
    }
    
    return 0;
}

static int my_log(PyObject *dict, Log *log) {
    assign_to_dict(dict, "score", log->score);
    assign_to_dict(dict, "episode_return", log->episode_return);
    assign_to_dict(dict, "episode_length", log->episode_length);
    assign_to_dict(dict, "wins", log->wins);
    assign_to_dict(dict, "losses", log->losses);
    assign_to_dict(dict, "draws", log->draws);
    assign_to_dict(dict, "timeouts", log->timeouts);
    assign_to_dict(dict, "invalid_action_rate", log->invalid_action_rate);
    assign_to_dict(dict, "chess_moves_completed", log->chess_moves_completed);
    float total = log->wins + log->losses + log->draws + log->timeouts;
    assign_to_dict(dict, "winrate", total > 0.0f ? log->wins / total : 0.0f);
    assign_to_dict(dict, "illegal_moves", log->illegal_moves);
    assign_to_dict(dict, "avg_legal_moves", log->avg_legal_moves);
    return 0;
}
