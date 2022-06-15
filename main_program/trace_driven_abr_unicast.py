import numpy as np
import scipy.io
import itertools
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--snr', type=str, help='low snr trace name')
args = parser.parse_args()

S_INFO = 5  # bit_rate, buffer_size, rebuffering_time, bandwidth_measurement, chunk_til_video_end
S_LEN = 8  # take how many frames in the past
MPC_FUTURE_CHUNK_COUNT = 1
M_IN_K = 1000.0
BUFFER_NORM_FACTOR = 10.0
CHUNK_TIL_VIDEO_END_CAP = 32.0
TOTAL_VIDEO_CHUNKS = 32
DEFAULT_QUALITY = 0  # default video quality without agent
REBUF_PENALTY = 4000  # 1 sec rebuffering -> this number of Mbps
SMOOTH_PENALTY = 0.5
TRAIN_SEQ_LEN = 100  # take as a train batch
MODEL_SAVE_INTERVAL = 100
RANDOM_SEED = 42
RAND_RANGE = 1000

CHUNK_COMBO_OPTIONS = []
# past errors in bandwidth
past_errors = []
past_bandwidth_ests = []

size_video1 = [1450544, 1456211, 1448768, 1475764, 1448377, 1450682, 1466778, 1444699, 1431871, 1432636, 1438689, 1398045, 1437565, 1347508, 1391073, 1387045, 1408341, 1382133, 1422908, 1401507, 1381447, 1393916, 1371175, 1391149, 1367623, 1364790, 1362530, 1376958, 1368421, 1387795, 1367348, 1372561, 1641125]
size_video2 = [523413, 532177, 529206, 533398, 525106, 529673, 526103, 522804, 518500, 521085, 520280, 508685, 515167, 471796, 493707, 481711, 507300, 495999, 510566, 497551, 502262, 506277, 501682, 506293, 502398, 499795, 497718, 502547, 501010, 508606, 503973, 503588, 571133]
size_video3 = [172298, 174137, 171884, 175931, 169597, 172825, 171121, 173919, 170134, 173417, 173125, 171710, 172919, 161674, 168840, 159595, 173027, 168332, 172122, 166160, 171588, 170332, 170704, 173129, 172148, 172648, 172451, 173907, 173569, 176446, 175209, 176819, 199137]
size_video4 = [81957, 80985, 79369, 80490, 80182, 80273, 79400, 82390, 80493, 83008, 81309, 80094, 83019, 79865, 82054, 79579, 83529, 80498, 82865, 80157, 83123, 82763, 81678, 83503, 84815, 84520, 85419, 85111, 84463, 85176, 85144, 85847, 96677]
VIDEO_BIT_RATE = [1000, 2000, 5000, 10000]  # kbps
bit_rate_map = {0: 1, 1: 2, 2: 5, 3: 10}
rest_time = 0.1-0.031*3
rest_time_mul = rest_time/3
ssim = []
ssim1 = []
ssim2 = []
ssim5 = []
ssim10 = []
f = open("ssim1.txt", "r")
for eachline in f.readlines():
    ssim1.append(float(eachline))
ssim.append(ssim1)
f = open("ssim2.txt", "r")
for eachline in f.readlines():
    ssim2.append(float(eachline))
ssim.append(ssim2)
f = open("ssim5.txt", "r")
for eachline in f.readlines():
    ssim5.append(float(eachline))
ssim.append(ssim5)
f = open("ssim10.txt", "r")
for eachline in f.readlines():
    ssim10.append(float(eachline))
ssim.append(ssim10)

psnr = []
psnr1 = []
psnr2 = []
psnr5 = []
psnr10 = []
f = open("psnr1.txt", "r")
for eachline in f.readlines():
    psnr1.append(float(eachline))
psnr.append(psnr1)
f = open("psnr2.txt", "r")
for eachline in f.readlines():
    psnr2.append(float(eachline))
psnr.append(psnr2)
f = open("psnr5.txt", "r")
for eachline in f.readlines():
    psnr5.append(float(eachline))
psnr.append(psnr5)
f = open("psnr10.txt", "r")
for eachline in f.readlines():
    psnr10.append(float(eachline))
psnr.append(psnr10)

s_batch = [np.zeros((S_INFO, S_LEN))]
last_bit_rate = DEFAULT_QUALITY
last_total_rebuf = 0
video_chunk_count = 0

input_dict = {'last_bit_rate': last_bit_rate,
            'last_total_rebuf': last_total_rebuf,
            'video_chunk_coount': video_chunk_count,
            's_batch': s_batch}
    
post_data = {'RebufferTime': 1000,
            'lastRequest': 0,
            'lastquality': 1,
            'lastthrp': 0.0,
            'buffer': 0.1}

# make chunk combination options
for combo in itertools.product([0,1,2,3], repeat=MPC_FUTURE_CHUNK_COUNT):
    CHUNK_COMBO_OPTIONS.append(combo)

def get_chunk_size(quality, index):
    if ( index < 0 or index > TOTAL_VIDEO_CHUNKS ):
        return 0
    # note that the quality and video labels are inverted (i.e., quality 8 is highest and this pertains to video1)
    sizes = {3: size_video1[index], 2: size_video2[index], 1: size_video3[index], 0: size_video4[index]}
    return sizes[quality]
    
def predict_future_bandwidth(updated_thrp):
    rebuffer_time = float(post_data['RebufferTime'] - input_dict['last_total_rebuf'])
    # retrieve previous state
    if len(input_dict['s_batch']) == 0:
        state = [np.zeros((S_INFO, S_LEN))]
    else:
        state = np.array(input_dict['s_batch'][-1], copy=True)

    # compute number of video chunks left
    video_chunk_remain = TOTAL_VIDEO_CHUNKS - input_dict['video_chunk_coount']
    input_dict['video_chunk_coount'] += 1

    # dequeue history record
    state = np.roll(state, -1, axis=1)

    # this should be S_INFO number of terms
    try:
        state[0, -1] = VIDEO_BIT_RATE[post_data['lastquality']] / float(np.max(VIDEO_BIT_RATE))
        state[1, -1] = post_data['buffer'] / BUFFER_NORM_FACTOR
        state[2, -1] = rebuffer_time / M_IN_K
        state[3, -1] = updated_thrp  # kilo byte / ms
        state[4, -1] = np.minimum(video_chunk_remain, CHUNK_TIL_VIDEO_END_CAP) / float(CHUNK_TIL_VIDEO_END_CAP)
        curr_error = 0 # defualt assumes that this is the first request so error is 0 since we have never predicted bandwidth
        if ( len(past_bandwidth_ests) > 0 ):
            curr_error  = abs(past_bandwidth_ests[-1]-state[3,-1])/float(state[3,-1])
        past_errors.append(curr_error)
    except ZeroDivisionError:
        # this should occur VERY rarely (1 out of 3000), should be a dash issue
        # in this case we ignore the observation and roll back to an eariler one
        past_errors.append(0)
        if len(input_dict['s_batch']) == 0:
            state = [np.zeros((S_INFO, S_LEN))]
        else:
            state = np.array(input_dict['s_batch'][-1], copy=True)


    # pick bitrate according to MPC           
    # first get harmonic mean of last 5 bandwidths
    past_bandwidths = state[3,-5:]
    while past_bandwidths[0] == 0.0:
        past_bandwidths = past_bandwidths[1:]
    #if ( len(state) < 5 ):
    #    past_bandwidths = state[3,-len(state):]
    #else:
    #    past_bandwidths = state[3,-5:]
    bandwidth_sum = 0
    for past_val in past_bandwidths:
        bandwidth_sum += (1/float(past_val))
    harmonic_bandwidth = 1.0/(bandwidth_sum/len(past_bandwidths))
    # future bandwidth prediction
    # divide by 1 + max of last 5 (or up to 5) errors
    max_error = 0
    error_pos = -5
    if ( len(past_errors) < 5 ):
        error_pos = -len(past_errors)
    max_error = float(max(past_errors[error_pos:]))
    future_bandwidth = harmonic_bandwidth/(1+max_error)
    past_bandwidth_ests.append(harmonic_bandwidth)
    input_dict['s_batch'].append(state)
    return future_bandwidth

def predict_future_bandwidth_fast(updated_thrp):
    rebuffer_time = float(post_data['RebufferTime'] - input_dict['last_total_rebuf'])
    # retrieve previous state
    if len(input_dict['s_batch']) == 0:
        state = [np.zeros((S_INFO, S_LEN))]
    else:
        state = np.array(input_dict['s_batch'][-1], copy=True)

    # compute number of video chunks left
    video_chunk_remain = TOTAL_VIDEO_CHUNKS - input_dict['video_chunk_coount']
    input_dict['video_chunk_coount'] += 1

    # dequeue history record
    state = np.roll(state, -1, axis=1)

    # this should be S_INFO number of terms
    try:
        state[0, -1] = VIDEO_BIT_RATE[post_data['lastquality']] / float(np.max(VIDEO_BIT_RATE))
        state[1, -1] = post_data['buffer'] / BUFFER_NORM_FACTOR
        state[2, -1] = rebuffer_time / M_IN_K
        state[3, -1] = updated_thrp  # kilo byte / ms
        state[4, -1] = np.minimum(video_chunk_remain, CHUNK_TIL_VIDEO_END_CAP) / float(CHUNK_TIL_VIDEO_END_CAP)
        curr_error = 0 # defualt assumes that this is the first request so error is 0 since we have never predicted bandwidth
        if ( len(past_bandwidth_ests) > 0 ):
            curr_error  = abs(past_bandwidth_ests[-1]-state[3,-1])/float(state[3,-1])
        past_errors.append(curr_error)
    except ZeroDivisionError:
        # this should occur VERY rarely (1 out of 3000), should be a dash issue
        # in this case we ignore the observation and roll back to an eariler one
        past_errors.append(0)
        if len(input_dict['s_batch']) == 0:
            state = [np.zeros((S_INFO, S_LEN))]
        else:
            state = np.array(input_dict['s_batch'][-1], copy=True)

    # pick bitrate according to MPC           
    # first get harmonic mean of last 5 bandwidths
    past_bandwidths = state[3,-5:]
    while past_bandwidths[0] == 0.0:
        past_bandwidths = past_bandwidths[1:]
    #if ( len(state) < 5 ):
    #    past_bandwidths = state[3,-len(state):]
    #else:
    #    past_bandwidths = state[3,-5:]
    bandwidth_sum = 0
    for past_val in past_bandwidths:
        bandwidth_sum += (1/float(past_val))
    future_bandwidth = 1.0/(bandwidth_sum/len(past_bandwidths))
    input_dict['s_batch'].append(state)
    return future_bandwidth

def run_abr(updated_thrp):
    # future_bandwidth = predict_future_bandwidth(updated_thrp)
    future_bandwidth = predict_future_bandwidth_fast(updated_thrp)
    print("future_bandwidth: ", future_bandwidth)
    # future chunks length (try 4 if that many remaining)
    last_index = int(post_data['lastRequest'])
    future_chunk_length = MPC_FUTURE_CHUNK_COUNT

    # all possible combinations of 5 chunk bitrates (9^5 options)
    # iterate over list and for each, compute reward and store max reward combination
    max_reward = -100000000
    best_combo = ()
    start_buffer = float(post_data['buffer'])
    #start = time.time()
    for full_combo in CHUNK_COMBO_OPTIONS:
        combo = full_combo[0:future_chunk_length]
        # calculate total rebuffer time for this combination (start with start_buffer and subtract
        # each download time and add 2 seconds in that order)
        curr_rebuffer_time = 0
        curr_buffer = rest_time
        bitrate_sum = 0
        smoothness_diffs = 0
        last_quality = int(post_data['lastquality'])
        for position in range(0, len(combo)):
            chunk_quality = combo[position]
            index = last_index + position + 1 # e.g., if last chunk is 3, then first iter is 3+0+1=4
            download_time = (get_chunk_size(chunk_quality, index)/1000000.)/(future_bandwidth/8.0) # this is MB/MB/s --> seconds
            if ( curr_buffer < download_time ):
                curr_rebuffer_time += (download_time - curr_buffer)
                curr_buffer = 0
            else:
                curr_buffer -= download_time
            curr_buffer += rest_time
            
            # linear reward
            bitrate_sum += VIDEO_BIT_RATE[chunk_quality]
            smoothness_diffs += abs(VIDEO_BIT_RATE[chunk_quality] - VIDEO_BIT_RATE[last_quality])

            last_quality = chunk_quality
        # compute reward for this combination (one reward per 5-chunk combo)
        # bitrates are in Mbits/s, rebuffer in seconds, and smoothness_diffs in Mbits/s
        
        # linear reward 
        reward = (bitrate_sum/1000.) - (REBUF_PENALTY*curr_rebuffer_time) # - (SMOOTH_PENALTY*smoothness_diffs/1000.)

        if ( reward > max_reward ):
            max_reward = reward
            best_combo = combo
    return best_combo[0]

# mat = scipy.io.loadmat("data_trace_driven_schedule_log_"+args.snr+".mat")
# mat = scipy.io.loadmat("data_trace_driven_schedule_log.mat")
mat = scipy.io.loadmat("data_trace_driven_schedule_log_env_move1.mat")
thrparr = mat["multicast"][:,4,2]

abr_ssim = []
abr_psnr = []
start = 0
quality = run_abr(thrparr[0])

avg_ssim = []
avg_psnr = []
last_thrp = thrparr[0]
for i in range(1, thrparr.shape[0]+1, 1):
    if i < thrparr.shape[0]:
        thrp = thrparr[i]
    else:
        thrp = last_thrp
    next_thrp = (thrp+last_thrp)/2.0  # average thrp between two 7 frames
    curr_quality = quality
    end = (start+1)%33
    size = get_chunk_size(curr_quality, start)
    download_time = float(size)/1000000/(last_thrp/8)
    eachssim = []
    eachpsnr = []
    if download_time < rest_time:
        if start == 32:
            eachssim.extend([ssim[curr_quality][start]]*4)
            eachpsnr.extend([psnr[curr_quality][start]]*4)
        else:
            eachssim.extend([ssim[curr_quality][start]]*3)
            eachpsnr.extend([psnr[curr_quality][start]]*3)
    else:
        sent_size = int((last_thrp/8)*rest_time*1000000)
        outname = "outframe"+str(start+1)+".yuv"
        size = str(sent_size)
        frame = str(start+1)
        rate = str(bit_rate_map[curr_quality])
        os.system("cp ../abr_videos2/frame%s_%s.h264 ./" % (frame, rate))
        os.system("truncate -s %s frame%s_%s.h264" % (size, frame, rate))
        os.system("ffmpeg -y -i frame%s_%s.h264 -c:v rawvideo -pix_fmt yuv420p outframe%s.yuv" % (frame, rate, frame))
        os.system("cp ../abr_videos2/frame%s_%s.h264 ./" % (frame, rate))
        file_size = os.path.getsize(outname)
        frame_size = int(4096*2160*1.5)
        if start == 32:
            gt_size = frame_size*4
        else:
            gt_size = frame_size*3
        diff_size = gt_size-file_size
        flag = False
        if (diff_size > 0):
            count = 0
            buffer = []
            with open(outname, "rb") as f:
                byte = f.read(frame_size)
                while byte != b"":
                    buffer.append(byte)
                    byte = f.read(frame_size)
            while diff_size > 0:
                # buffer.append(buffer[0])  # reuse the last frame if the current frame is lost
                buffer.append(bytearray(frame_size))  # fill it with a blank frame
                diff_size -= frame_size
                flag = True
            f = open(outname, 'wb')
            for eachbuf in buffer:
                f.write(eachbuf)
            f.close()

        os.system("ffmpeg -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i ../abr_videos2/frame%s.yuv -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i outframe%s.yuv -lavfi psnr=stats_file=log.txt -f null -" % (frame, frame))
        os.system("grep -oP '(?<=psnr_avg:)[0-9]+.[0-9]+' log.txt > psnr.txt") 
        os.system("ffmpeg -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i ../abr_videos2/frame%s.yuv -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i outframe%s.yuv -lavfi ssim=stats_file=log.txt -f null -" % (frame, frame))
        os.system("grep -oP '(?<=All:)[0-9].[0-9]+' log.txt > ssim.txt")
        f = open("psnr.txt")
        for eachline in f.readlines():
            eachpsnr.append(float(eachline))
        f.close()
        f = open("ssim.txt")
        for eachline in f.readlines():
            eachssim.append(float(eachline))
        f.close()

    start = end
    end = (start+1)%33
    size = get_chunk_size(curr_quality, start)
    download_time = float(size)/1000000/(next_thrp/8)
    if download_time < rest_time:
        if start == 32:
            eachssim.extend([ssim[curr_quality][start]]*4)
            eachpsnr.extend([psnr[curr_quality][start]]*4)
        else:
            eachssim.extend([ssim[curr_quality][start]]*3)
            eachpsnr.extend([psnr[curr_quality][start]]*3)
    else:
        sent_size = int((next_thrp/8)*rest_time*1000000)
        outname = "outframe"+str(start+1)+".yuv"
        size = str(sent_size)
        frame = str(start+1)
        rate = str(bit_rate_map[curr_quality])
        os.system("cp ../abr_videos2/frame%s_%s.h264 ./" % (frame, rate))
        os.system("truncate -s %s frame%s_%s.h264" % (size, frame, rate))
        os.system("ffmpeg -y -i frame%s_%s.h264 -c:v rawvideo -pix_fmt yuv420p outframe%s.yuv" % (frame, rate, frame))
        os.system("cp ../abr_videos2/frame%s_%s.h264 ./" % (frame, rate))
        file_size = os.path.getsize(outname)
        frame_size = int(4096*2160*1.5)
        if start == 32:
            gt_size = frame_size*4
        else:
            gt_size = frame_size*3
        diff_size = gt_size-file_size
        if (diff_size > 0):
            f = open(outname, 'rb')
            count = 0
            buffer = []
            with open(outname, "rb") as f:
                byte = f.read(frame_size)
                while byte != b"":
                    buffer.append(byte)
                    byte = f.read(frame_size)
            while diff_size > 0:
                # buffer.append(buffer[0])  # reuse the last frame if the current frame is lost
                buffer.append(bytearray(frame_size))  # fill it with a blank frame
                diff_size -= frame_size
            f = open(outname, 'wb')
            for eachbuf in buffer:
                f.write(eachbuf)
            f.close()
        os.system("ffmpeg -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i ../abr_videos2/frame%s.yuv -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i outframe%s.yuv -lavfi psnr=stats_file=log.txt -f null -" % (frame, frame))
        os.system("grep -oP '(?<=psnr_avg:)[0-9]+.[0-9]+' log.txt > psnr.txt")
        os.system("ffmpeg -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i ../abr_videos2/frame%s.yuv -f rawvideo -s 4096x2160 -pix_fmt yuv420p -i outframe%s.yuv -lavfi ssim=stats_file=log.txt -f null -" % (frame, frame))
        os.system("grep -oP '(?<=All:)[0-9].[0-9]+' log.txt > ssim.txt")
        f = open("psnr.txt")
        for eachline in f.readlines():
            eachpsnr.append(float(eachline))
        f.close()
        f = open("ssim.txt")
        for eachline in f.readlines():
            eachssim.append(float(eachline))
        f.close()
 
    while len(eachpsnr) < 7:
        eachpsnr.append(eachpsnr[-1])
    avg_psnr.extend(eachpsnr)
    while len(eachssim) < 7:
        eachssim.append(eachssim[-1])
    avg_ssim.extend(eachssim)
    input_dict['last_bit_rate'] = VIDEO_BIT_RATE[curr_quality]
    post_data['lastRequest'] += 1
    post_data['lastRequest'] = post_data['lastRequest'] % 33
    post_data['lastquality'] = curr_quality
    post_data['lastthrp'] = last_thrp
    print("i: ", i, ", curr_quality: ", curr_quality)
    quality = run_abr(last_thrp)
    
    last_thrp = thrp

abr_ssim = avg_ssim
abr_psnr = avg_psnr

# with open("fast_ssim_"+args.snr+"_unicast"+".txt", "w") as output:
# with open("fast_ssim_old_unicast_reuse.txt", "w") as output:
with open("fast_ssim_env_unicast"+".txt", "w") as output:
    for i in range(len(abr_ssim)):
        output.write(str(abr_ssim[i]))
        output.write("\n")

# with open("fast_psnr_"+args.snr+"_unicast"+".txt", "w") as output:
# with open("fast_psnr_old_unicast_reuse.txt", "w") as output:
with open("fast_psnr_env_unicast"+".txt", "w") as output:
    for i in range(len(abr_psnr)):
        output.write(str(abr_psnr[i]))
        output.write("\n")

# print("final_ssim: ", abr_ssim)
# print("final_psnr: ", abr_psnr)
