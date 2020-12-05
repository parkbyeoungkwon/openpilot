import os
import fcntl
import signal
import json
import time

LIMIT_PATH = '/data/data/com.neokii.oproadlimit/files/'
LIMIT_FILE = '/data/data/com.neokii.oproadlimit/files/oproadlimit.json'

current_milli_time = lambda: int(round(time.time() * 1000))


class RoadSpeedLimiter:
  def __init__(self):
    self.json = None
    self.last_updated = 0

    try:
      os.remove(LIMIT_FILE)
    except:
      pass

    try:
      signal.signal(signal.SIGIO, self.handler)
      fd = os.open(LIMIT_PATH, os.O_RDONLY)
      fcntl.fcntl(fd, fcntl.F_SETSIG, 0)
      fcntl.fcntl(fd, fcntl.F_NOTIFY, fcntl.DN_MODIFY | fcntl.DN_CREATE | fcntl.DN_MULTISHOT)
    except Exception as ex:
      pass

  def handler(self, signum, frame):
    try:
      self.json = None
      if os.path.isfile(LIMIT_FILE):
        with open(LIMIT_FILE, 'r') as f:
          self.json = json.load(f)
          self.last_updated = current_milli_time()

    except Exception as ex:
      pass

  def get_val(self, key):
    if key in self.json:
      return self.json[key]
    return None

  def get_max_speed(self, CS, v_cruise_kph):

    if current_milli_time() - self.last_updated > 1000 * 10:
      return 0, 0

    try:

      # car_speed_kph = CS.vEgo * 3.6

      road_limit_speed = self.get_val('road_limit_speed')
      is_highway = self.get_val('is_highway')

      cam_limit_speed_left_dist = self.get_val('cam_limit_speed_left_dist')
      cam_limit_speed = self.get_val('cam_limit_speed')

      section_limit_speed = self.get_val('section_limit_speed')
      # section_avg_speed = self.get_val('section_avg_speed')
      section_left_dist = self.get_val('section_left_dist')
      # section_left_time = self.get_val('section_left_time')

      if is_highway is not None:
        if is_highway:
          MIN_LIMIT = 60
          MAX_LIMIT = 120
        else:
          MIN_LIMIT = 30
          MAX_LIMIT = 100
      else:
        MIN_LIMIT = 30
        MAX_LIMIT = 120

      if cam_limit_speed_left_dist is not None and cam_limit_speed is not None and cam_limit_speed_left_dist > 0:
        if MIN_LIMIT <= cam_limit_speed <= MAX_LIMIT and cam_limit_speed_left_dist < (cam_limit_speed / 3.6) * 10:
          return cam_limit_speed, cam_limit_speed

        return 0, cam_limit_speed

      elif section_left_dist is not None and section_limit_speed is not None and section_left_dist > 0:
        if MIN_LIMIT <= section_limit_speed <= MAX_LIMIT:
          return section_limit_speed, section_limit_speed

        return 0, section_limit_speed

    except:
      pass

    return 0, 0


road_speed_limiter = None


def road_speed_limiter_get_max_speed(CS, v_cruise_kph):
  global road_speed_limiter
  if road_speed_limiter is None:
    road_speed_limiter = RoadSpeedLimiter()

  return road_speed_limiter.get_max_speed(CS, v_cruise_kph)
