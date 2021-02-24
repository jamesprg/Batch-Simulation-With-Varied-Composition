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
        for j in temp:
            #if the h-bond matches the current set
            if ((i == j[1:4]).all(-1).any()):
                #add the interaction to this combos array
                frames.append(j[0])
        #add frames to this combos interaction
        if n == 0:
            interactions = []
            interactions.append(frames)
        else:
            interactions.append(frames)
        # we could skip the interactions list
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
            #if the current frame minus the last saved frame is < 4, cont
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