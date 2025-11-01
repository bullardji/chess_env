import gymnasium
import numpy as np

import pufferlib
from pufferlib.ocean.chess import binding

class Chess(pufferlib.PufferEnv):
    def __init__(self, num_envs=1, render_mode=None, log_interval=128, buf=None, seed=0,
                 max_moves=3000, opponent_depth=-1, reward_shaping_weight=0.05,
                 reward_draw=-0.5, reward_invalid_piece=-0.1, reward_invalid_move=-0.1, 
                 reward_valid_piece=0.0, reward_valid_move=0.0, render_fps=30, human_play=0,
                 starting_fen="rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1"):
        self.single_observation_space = gymnasium.spaces.Box(
            low=0, high=255, shape=(1045,), dtype=np.uint8)
        self.single_action_space = gymnasium.spaces.Discrete(64)
        self.render_mode = render_mode
        self.num_agents = num_envs
        self.log_interval = log_interval
        self.tick = 0
        
        super().__init__(buf)
        self.c_envs = binding.vec_init(
            self.observations, self.actions, self.rewards,
            self.terminals, self.truncations, num_envs, seed,
            max_moves=max_moves, opponent_depth=opponent_depth,
            reward_shaping_weight=reward_shaping_weight, reward_draw=reward_draw,
            reward_invalid_piece=reward_invalid_piece, reward_invalid_move=reward_invalid_move,
            reward_valid_piece=reward_valid_piece, reward_valid_move=reward_valid_move,
            render_fps=render_fps, human_play=human_play, starting_fen=starting_fen)
    
    def reset(self, seed=0):
        binding.vec_reset(self.c_envs, seed)
        self.tick = 0
        return self.observations, []
    
    def step(self, actions):
        self.tick += 1
        self.actions[:] = actions
        binding.vec_step(self.c_envs)
        info = [binding.vec_log(self.c_envs)] if self.tick % self.log_interval == 0 else []
        return self.observations, self.rewards, self.terminals, self.truncations, info
    
    def render(self):
        binding.vec_render(self.c_envs, 0)
    
    def close(self):
        binding.vec_close(self.c_envs)

if __name__ == '__main__':
    N = 4096
    env = Chess(num_envs=N)
    env.reset()
    steps = 0

    CACHE = 1024
    actions = np.random.randint(0, 64, (CACHE, N))

    import time
    start = time.time()
    while time.time() - start < 10:
        env.step(actions[steps % CACHE])
        steps += 1

    print('Chess SPS:', int(env.num_agents * steps / (time.time() - start)))
