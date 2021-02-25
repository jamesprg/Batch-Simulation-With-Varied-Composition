import numpy as np

def hbond_lifetimes(file,intermittency=0,singles=False):
    temp_times = []
    #get the set combinations
    for n,i in enumerate(file):
        if n == 0:
            combos = i[1:4]
        elif ((combos == i[1:4]).all(-1).any()) == False:
            combos = np.vstack((combos,i[1:4]))
    #now we want to make an array for each id set
    for n,i in enumerate(combos):
        frames = []
        for j in file:
            #if the h-bond matches the current set
            if ((i == j[1:4]).all(-1).any()):
                #add the interaction to this combos array
                frames.append(j[0])
        #get the life time of each interaction
        start = 0
        for fn, f in enumerate(frames):
            if f == frames[0]:
                start = f
            #check for end of frames
            elif frames[fn] == frames[-1]:
                time = frames[fn-1]-start
                if singles:
                    temp_times.append(time)
                else:
                    if time != 0 :
                        temp_times.append(time)
                continue
            #if on the next frame (or within intermittency)
            elif (f-frames[fn-1]) <= (1+intermittency):
                continue
            else:
                time = frames[fn-1]-start
                if singles:
                    temp_times.append(time)
                else:
                    if time != 0:
                        temp_times.append(time)
                #update start of interaction
                start = f
    return temp_times

def fast_viterbi(obs, states, start_p, trans_p, emission_p):
    """Algorithm used to decode a signal and find the optimal path.

    Parameters
    ----------
    obs : numpy.ndarray
        Numerical representations of the observations.
    states : list of int
        Numerical representations of the states in the system.
    start_p : list of float
        Probabilities of starting in a particular hidden state.
    trans_p : list of list of float
        Probabilities of transitioning from one hidden state to
        another.
    emission_p : list of list of float
        Probabilities of emitting an observable given the present
        hidden state.

    Returns
    -------
    optimal_path : list of int
        Optimal path of the system given the observations and Hidden
        Markov Model parameters."""

    V = np.zeros((len(obs), len(states)))
    prev_states = np.zeros((len(obs), len(states)), dtype=int)

    V[0, :] = np.log10(start_p[:] * emission_p[:, obs[0]])

    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        for st in states:
            max_tr_prob = (V[t - 1, :] + np.log10(trans_p[:, st])).max()
            V[t, st] = max_tr_prob + np.log10(emission_p[st, obs[t]])
            prev_states[t, st] = np.where(V[t - 1, :] +
                                          np.log10(trans_p[:, st]) ==
                                          max_tr_prob)[0][0]

    optimal_path = np.zeros(len(obs), dtype=int)

    # The highest final probability.
    max_final_prob = V[-1, :].max()

    # Get most probable state and its backtrack.
    optimal_path[-1] = np.where(V[-1, :] == max_final_prob)[0][0]
    previous = optimal_path[-1]

    # Follow the each backtrack from the saved paths.
    for t in range(len(V) - 2, -1, -1):
        optimal_path[t] = prev_states[t + 1, previous]
        previous = optimal_path[t]

    return optimal_path

def hbond_signal(file, states = [0, 1],starting = [0.9, 0.1],
                transition = [[0.7, 0.3],[0.3, 0.7]], 
                emission = [[0.7, 0.3],[0.3, 0.7]]):
    states = np.array(states)
    start_p = np.array(starting)
    trans_p = np.array(transition)
    emission_p = np.array(emission)
    temp_times = []
    #get the set combinations
    for n,i in enumerate(file):
        if n == 0:
            combos = i[1:4]
        elif ((combos == i[1:4]).all(-1).any()) == False:
            combos = np.vstack((combos,i[1:4]))
    #now we want to make an array for each id set
    for n,i in enumerate(combos):
        frames = []
        for j in file:
            #if the h-bond matches the current set
            if ((i == j[1:4]).all(-1).any()):
                #add the interaction to this combos array
                frames.append(j[0])
        #clean frames using viterbi algo.
        joined = np.zeros(int(frames[-1]+1))
        for i in frames:
            joined[int(i)] = 1
        joined = joined.astype(int)
        cleaned = fast_viterbi(joined, states, start_p, trans_p, emission_p)
        frames = np.where(cleaned==1)[0]
        #get the life time of each interaction
        start = 0
        for fn, f in enumerate(frames):
            if f == frames[0]:
                start = f
            #check for end of frames
            elif frames[fn] == frames[-1]:
                time = frames[fn-1]-start
                if time != 0 :
                        temp_times.append(time)
                continue
            #if on the next frame (or within intermittency)
            elif (f-frames[fn-1]) <= 1:
                continue
            else:
                time = frames[fn-1]-start
                if time != 0:
                        temp_times.append(time)
                #update start of interaction
                start = f
    return temp_times

