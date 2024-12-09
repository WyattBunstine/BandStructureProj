import os

import numpy as np
import cv2


def find_images(path):
    images = []
    for file in os.listdir(path):
        if os.path.isdir(path + file):
            images = images + find_images(path + file + "/")
        elif file.endswith(".png"):
            if ("data" in path and "venv" not in path and ".git" not in path and ".idea" not in path and "candidates" in path):
                images.append(path + file)
    return images


def make_collage():
    file_loc = "/home/wyatt/PycharmProjects/"
    image_size = 500
    images = find_images(file_loc)
    print(len(images))
    size = int(np.sqrt(len(images)))
    final_image = None
    for i in np.arange(size):
        horiz_image = cv2.imread(images[size * i])
        horiz_image = cv2.resize(horiz_image, (image_size*2, image_size))
        for j in np.arange(1, size):
            if size * i + j < len(images):
                im = cv2.imread(images[size * i + j])
                im = cv2.resize(im, (image_size*2, image_size))
                horiz_image = np.hstack((horiz_image, im))
            else:
                horiz_image = np.hstack((horiz_image, np.zeros((image_size, image_size*2, 3))))
        if final_image is None:
            final_image = horiz_image
        else:
            final_image = np.vstack((final_image, horiz_image))
    final_image = cv2.resize(final_image, (1024*2, 1024))
    cv2.imwrite("candidate_collage.png", final_image)


if __name__ == '__main__':
    make_collage()
