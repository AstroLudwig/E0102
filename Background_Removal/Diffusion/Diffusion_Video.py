import cv2
import os



def create_video(image_prefix_filename,save_filename,start,stop):
	prefix = image_prefix_filename
	video_name = save_filename
	images = []
	for i in range(start,stop,5):
		#images.append(img_pref+str(i)+img_suff)
		images.append(prefix+str(i)+'.png')

	#frame = cv2.imread(os.path.join(image_folder,images[0]))
	frame = cv2.imread(images[0])
	height,width,layers=frame.shape

	fourcc = cv2.VideoWriter_fourcc(*'MJPG')

	video = cv2.VideoWriter(video_name,fourcc,40,(width,height))
	count = 5
	for image in images:
		print('Writing img: '+str(count))
		count += 5
		#video.write(cv2.imread(os.path.join(image_folder,image)))
		video.write(cv2.imread(image))

	cv2.destroyAllWindows()
	video.release()	

"""
baseLength = 9
targetLength = baseLength+4
filename = "im24_25.png"

def insert0(position, file) :
	return file[:position] + '0' + file[position:]

image_folder = 'im24/imgs/'

images = [img for img in os.listdir(image_folder)]

print(len(images))
print()

print(insert0(5,filename))
"""