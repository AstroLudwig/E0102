# -*- Copyright (c) 2018, Bethany Ann Ludwig, All rights reserved. -*-
"""
NAME:
    Diffusion Video
PURPOSE:
    Module for creating a video from png images.
"""

import cv2
import os

def create_video(image_prefix_filename,save_filename,start,stop):
	prefix = image_prefix_filename
	video_name = save_filename
	images = []
	for i in range(start,stop,5):
		images.append(prefix+str(i)+'.png')

	frame = cv2.imread(images[0])
	height,width,layers=frame.shape

	fourcc = cv2.VideoWriter_fourcc(*'MJPG')

	video = cv2.VideoWriter(video_name,fourcc,40,(width,height))
	count = 5
	for image in images:
		print('Writing img: '+str(count))
		count += 5
		video.write(cv2.imread(image))

	cv2.destroyAllWindows()
	video.release()	
